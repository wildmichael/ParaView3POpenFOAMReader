/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkOpenFOAMReader.cxx,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// Thanks to Terry Jordan of SAIC at the National Energy
// Technology Laboratory who developed this class.
// Please address all comments to Terry Jordan (terry.jordan@sa.netl.doe.gov)
//
// Token-based FoamFile format lexer/parser, performance enhancements,
// gzipped file and lagrangian field support by Takuya Oshima
// (oshima@eng.niigata-u.ac.jp)
//
// * GUI Based selection of mesh regions and fields available in the case
// * Minor bug fixes / Strict memory allocation checks
// * Minor performance enhancements
// by Philippose Rajan (sarith@rocketmail.com)

// version 2008-08-30

// Hijack the CRC routine of zlib to omit CRC check for gzipped files
// (on OSes other than Windows where the mechanism doesn't work due
// to pre-bound DLL symbols) if set to 1, or not (set to 0). Affects
// performance by about 3% - 4%.
#define VTK_FOAMFILE_OMIT_CRCCHECK 1

// The input/output buffer sizes for zlib in bytes.
#define VTK_FOAMFILE_INBUFSIZE (16384)
#define VTK_FOAMFILE_OUTBUFSIZE (131072)
#define VTK_FOAMFILE_INCLUDE_STACK_SIZE (10)

#if defined(_MSC_VER) && (_MSC_VER >= 1400)
#define _CRT_SECURE_NO_WARNINGS 1
#endif

#if VTK_FOAMFILE_OMIT_CRCCHECK
#define ZLIB_INTERNAL
#endif

// for possible future extension of linehead-aware directives
#define VTK_FOAMFILE_RECOGNIZE_LINEHEAD 0

#include "vtkOpenFOAMReader.h"

#include <vtkstd/vector>
#include "vtksys/DateStamp.h"
#include "vtksys/SystemTools.hxx"
#include <vtksys/ios/sstream>
#include <vtkzlib/zlib.h>

#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkCharArray.h"
#include "vtkCollection.h"
#include "vtkConvexPointSet.h"
#include "vtkDataArraySelection.h"
#include "vtkDirectory.h"
#include "vtkDoubleArray.h"
#include "vtkFloatArray.h"
#include "vtkHexahedron.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkIntArray.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkPolygon.h"
#include "vtkPyramid.h"
#include "vtkQuad.h"
#include "vtkSortDataArray.h"
#include "vtkStdString.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkStringArray.h"
#include "vtkTetra.h"
#include "vtkTriangle.h"
#include "vtkUnstructuredGrid.h"
#include "vtkVertex.h"
#include "vtkWedge.h"

#if !(defined(_WIN32) && !defined(__CYGWIN__) || defined(__LIBCATAMOUNT__))
// for getpwnam() / getpwuid()
#include <sys/types.h>
#include <pwd.h>
// for getuid()
#include <unistd.h>
#endif
// for fabs()
#include <math.h>

#if VTK_FOAMFILE_OMIT_CRCCHECK
uLong ZEXPORT crc32(uLong, const Bytef *, uInt)
{ return 0; }
#endif

vtkCxxRevisionMacro(vtkOpenFOAMReader, "$Revision: 1.13 $");
vtkStandardNewMacro(vtkOpenFOAMReader);

//-----------------------------------------------------------------------------
// class vtkOpenFOAMReaderPrivate
// the reader core of vtkOpenFOAMReader
class VTK_IO_EXPORT vtkOpenFOAMReaderPrivate : public vtkObject
{
public:
  static vtkOpenFOAMReaderPrivate *New();
  vtkTypeRevisionMacro(vtkOpenFOAMReaderPrivate, vtkObject);

  vtkDoubleArray *GetTimeValues() { return this->TimeValues; }
  vtkGetMacro(TimeStep, int);
  vtkSetMacro(TimeStep, int);
  const vtkStdString &GetRegionName() const { return this->RegionName; }

  // gather timestep information
  bool MakeInformationVector(const vtkStdString &, const vtkStdString &,
    const vtkStdString &, vtkOpenFOAMReader *);
  // read mesh/fields and create dataset
  int RequestData(vtkMultiBlockDataSet *, bool, bool, bool);
  void SetTimeValue(const double);
  int MakeMetaDataAtTimeStep(vtkStringArray *, vtkStringArray *,
    vtkStringArray *, const bool);
  void SetupInformation(const vtkStdString &, const vtkStdString &,
    const vtkStdString &, vtkOpenFOAMReaderPrivate *);

private:
  struct intArrayVector;
  struct intVectorVector;

  struct vtkFoamError;
  struct vtkFoamToken;
  struct vtkFoamFileStack;
  struct vtkFoamFile;
  struct vtkFoamIOobject;
  struct vtkFoamEntryValue;
  struct vtkFoamEntry;
  struct vtkFoamDict;

  struct vtkFoamBoundaryEntry
  {
    enum bt
      {
      PHYSICAL = 1, // patch, wall
      PROCESSOR = 2, // processor
      GEOMETRICAL = 0, // symmetryPlane, wedge, cyclic, empty, etc.
      };
    vtkStdString boundaryName;
    int nFaces, startFace, allBoundariesStartFace;
    bool isActive;
    bt boundaryType;
  };

  struct vtkFoamBoundaryDict : public vtkstd::vector<vtkFoamBoundaryEntry>
  {
    // we need to keep the path to time directory where the current mesh
    // is read from, since boundaryDict may be accessed multiple times
    // at a timestep for patch selections
    vtkStdString timeDir;
  };

  vtkOpenFOAMReader *Parent;

  // case and region
  vtkStdString CasePath;
  vtkStdString RegionName;
  vtkStdString ProcessorName;

  // time information
  vtkDoubleArray *TimeValues;
  int TimeStep;
  int TimeStepOld;
  vtkStringArray *TimeNames;

  int InternalMeshSelectionStatus;
  int InternalMeshSelectionStatusOld;

  // filenames / directories
  vtkStringArray *VolFieldFiles;
  vtkStringArray *PointFieldFiles;
  vtkStringArray *LagrangianFieldFiles;
  vtkStringArray *PolyMeshPointsDir;
  vtkStringArray *PolyMeshFacesDir;

  // for mesh construction
  vtkIdType NumCells;
  vtkIdType NumPoints;
  vtkIntArray *FaceOwner;

  // for cell-to-point interpolation
  vtkUnstructuredGrid *AllBoundaries;
  vtkIntArray *AllBoundariesPointMap;
  vtkIntArray *InternalPoints;

  // for caching mesh
  vtkUnstructuredGrid *InternalMesh;
  vtkMultiBlockDataSet *BoundaryMesh;
  intArrayVector *BoundaryPointMap;
  vtkFoamBoundaryDict BoundaryDict;
  vtkMultiBlockDataSet *PointZoneMesh;
  vtkMultiBlockDataSet *FaceZoneMesh;
  vtkMultiBlockDataSet *CellZoneMesh;

  // for polyhedra handling
  vtkIntArray *AdditionalCellIds;
  intArrayVector *AdditionalCellPoints;

  // constructor and destructor are kept private
  vtkOpenFOAMReaderPrivate();
  ~vtkOpenFOAMReaderPrivate();

  // not implemented.
  vtkOpenFOAMReaderPrivate(const vtkOpenFOAMReaderPrivate &);
  void operator=(const vtkOpenFOAMReaderPrivate &);

  // clear mesh construction
  void ClearInternalMeshes();
  void ClearBoundaryMeshes();
  void ClearMeshes();

  vtkStdString RegionPath() const
  { return (this->RegionName == "" ? "" : "/")  + this->RegionName; }
  vtkStdString TimePath(const int timeI) const
  { return this->CasePath + this->TimeNames->GetValue(timeI); }
  vtkStdString TimeRegionPath(const int timeI) const
  { return this->TimePath(timeI) + this->RegionPath(); }
  vtkStdString CurrentTimePath() const
  { return this->TimePath(this->TimeStep); }
  vtkStdString CurrentTimeRegionPath() const
  { return this->TimeRegionPath(this->TimeStep); }
  vtkStdString CurrentTimeRegionMeshPath(vtkStringArray *dir) const
  { return this->CasePath + dir->GetValue(this->TimeStep) + this->RegionPath()
    + "/polyMesh/"; }
  vtkStdString RegionPrefix() const
  { return this->RegionName + (this->RegionName == "" ? "" : "/"); }

  // search time directories for mesh
  void AppendMeshDirToArray(vtkStringArray *, const vtkStdString &, const int);
  void PopulatePolyMeshDirArrays();

  // search a time directory for field objects
  void GetFieldNames(const vtkStdString &, const bool, vtkStringArray *,
    vtkStringArray *);
  void SortFieldFiles(vtkStringArray *, vtkStringArray *, vtkStringArray *);
  void LocateLagrangianClouds(vtkStringArray *, const vtkStdString &);

  // read controlDict
  bool ListTimeDirectoriesByControlDict(vtkFoamDict *dict);
  bool ListTimeDirectoriesByInstances();

  // read mesh files
  vtkFloatArray* ReadPointsFile();
  intVectorVector* ReadFacesFile (const vtkStdString &);
  intVectorVector* ReadOwnerNeighborFiles(const vtkStdString &,
    intVectorVector *);
  bool CheckFacePoints(intVectorVector *);

  // create mesh
  void InsertCellsToGrid(vtkUnstructuredGrid *, const intVectorVector *,
    const intVectorVector *, vtkFloatArray *, vtkIdTypeArray *,
    vtkIntArray *);
  vtkUnstructuredGrid *MakeInternalMesh(const intVectorVector *,
    const intVectorVector *, vtkFloatArray *);
  void InsertFacesToGrid(vtkUnstructuredGrid *, const intVectorVector *,
    int, int, vtkIntArray *, vtkIdList *, vtkIntArray *, const bool);
  vtkMultiBlockDataSet* MakeBoundaryMesh(const intVectorVector *,
    vtkFloatArray *);
  void SetBlockName(vtkMultiBlockDataSet *, unsigned int, const char *);
  void TruncateFaceOwner();

  // move additional points for decomposed cells
  vtkPoints *MoveInternalMesh(vtkUnstructuredGrid *, vtkFloatArray *);
  void MoveBoundaryMesh(vtkMultiBlockDataSet *, vtkFloatArray *);

  // cell-to-point interpolator
  void InterpolateCellToPoint(vtkFloatArray *, vtkFloatArray *,
    vtkUnstructuredGrid *, vtkIntArray *, const int);

  // read and create cell/point fields
  void ConstructDimensions(vtkStdString *, vtkFoamDict *);
  bool ReadFieldFile(vtkFoamIOobject *, vtkFoamDict *, const vtkStdString &,
    vtkDataArraySelection *);
  vtkFloatArray *FillField(vtkFoamEntry *, int, vtkFoamIOobject *,
    const vtkStdString &);
  void GetVolFieldAtTimeStep(vtkUnstructuredGrid *, vtkMultiBlockDataSet *,
    const vtkStdString &);
  void GetPointFieldAtTimeStep(vtkUnstructuredGrid *, vtkMultiBlockDataSet *,
    const vtkStdString &);
  void AddArrayToFieldData(vtkDataSetAttributes *, vtkDataArray *,
    const vtkStdString &);

  // create lagrangian mesh/fields
  vtkMultiBlockDataSet *MakeLagrangianMesh();

  // create point/face/cell zones
  vtkFoamDict *GatherBlocks(const char *, bool);
  bool GetPointZoneMesh(vtkMultiBlockDataSet *, vtkPoints *);
  bool GetFaceZoneMesh(vtkMultiBlockDataSet *, const intVectorVector *,
    vtkPoints *);
  bool GetCellZoneMesh(vtkMultiBlockDataSet *, const intVectorVector *,
    const intVectorVector *, vtkPoints *);
};

vtkCxxRevisionMacro(vtkOpenFOAMReaderPrivate, "$Revision: 1.00 $");
vtkStandardNewMacro(vtkOpenFOAMReaderPrivate);

//-----------------------------------------------------------------------------
// struct intArrayVector
struct vtkOpenFOAMReaderPrivate::intArrayVector
  : public vtkstd::vector<vtkIntArray *>
{
private:
  typedef vtkstd::vector<vtkIntArray *> Superclass;

public:
  ~intArrayVector()
  {
    for(size_t arrayI = 0; arrayI < Superclass::size(); arrayI++)
      {
      if(Superclass::operator[](arrayI) != NULL)
        {
        Superclass::operator[](arrayI)->Delete();
        }
      }
  }
};

//-----------------------------------------------------------------------------
// struct intVectorVector
struct vtkOpenFOAMReaderPrivate::intVectorVector
{
private:
  vtkIntArray *body_, *indices_;
  int nElements_;
  void clear() { indices_->Delete(); body_->Delete(); }
  intVectorVector();

public:
  intVectorVector(const intVectorVector &ivv)
    : body_(ivv.body_), indices_(ivv.indices_),
      nElements_(ivv.nElements_)
  {
    indices_->Register(0); // vtkDataArrays do not have ShallowCopy
    body_->Register(0);
  }
  intVectorVector(const int nElements, const int bodyLength)
    : body_(vtkIntArray::New()), indices_(vtkIntArray::New()),
      nElements_(nElements)
  {
    indices_->SetNumberOfValues(nElements + 1);
    body_->SetNumberOfValues(bodyLength);
  }

  ~intVectorVector() { clear(); }

  // GetSize() returns all allocated size (Size) while GetDataSize()
  // returns actually used size (MaxId * nComponents)
  int bodySize() const { return body_->GetSize(); }
  // note that vtkIntArray::Resize() allocates (current size + new
  // size) bytes if current size < new size
  void resizeBody(const int bodyLength) { body_->Resize(bodyLength); }
  int *setIndex(const int i, const int bodyI)
  { return body_->GetPointer(*indices_->GetPointer(i) = bodyI); }
  void setValue(const int bodyI, int value) { body_->SetValue(bodyI, value); }
  const int *operator[](const int i) const
  { return body_->GetPointer(indices_->GetValue(i)); }
  int size(const int i) const
  { return indices_->GetValue(i + 1) - indices_->GetValue(i); }
  int nElements() const { return nElements_; }
  vtkIntArray *indices() { return indices_; }
  vtkIntArray *body() { return body_; }
};

//-----------------------------------------------------------------------------
// class vtkFoamError
// class for exception-carrying object
struct vtkOpenFOAMReaderPrivate::vtkFoamError: public vtkStdString
{
public:
  vtkFoamError(): vtkStdString() {}
  vtkFoamError(const vtkFoamError& e): vtkStdString(e) {}
  ~vtkFoamError() {}
  // a super-easy way to make use of operator<<()'s defined in
  // vtksys_ios::ostringstream class
  template <class T> vtkFoamError& operator<<(const T& t)
  { vtksys_ios::ostringstream os; os << t; operator+=(os.str()); return *this; }
};

//-----------------------------------------------------------------------------
// class vtkFoamToken
// token class which also works as container for list types
// - a word token is treated as a string token for simplicity
// - handles only atomic types. Handling of list types are left to the
//   derived classes.
struct vtkOpenFOAMReaderPrivate::vtkFoamToken
{
public:
  enum tokenType
  {
    // undefined type
    UNDEFINED,
    // atomic types
    PUNCTUATION, LABEL, SCALAR, STRING, IDENTIFIER,
    // vtkObject-derived list types
    STRINGLIST, LABELLIST, SCALARLIST, VECTORLIST,
    // original list types
    LABELLISTLIST, ENTRYVALUELIST, EMPTYLIST, DICTIONARY,
    // exceptional state
    UNIFORMLABELLIST, UNIFORMSCALARLIST,
    // error state
    ERROR
  };

protected:
  tokenType type_;
  union
  {
    char char_;
    int int_;
    double double_;
    vtkStdString* string_;
    vtkObjectBase *vtkObjectPtr_;
    // vtkObject-derived list types
    vtkIntArray *labelListPtr_;
    vtkFloatArray *scalarListPtr_, *vectorListPtr_;
    vtkStringArray *stringListPtr_;
    // original list types
    intVectorVector *labelListListPtr_;
    vtkstd::vector<vtkFoamEntryValue*> *entryValuePtrs_;
    vtkFoamDict *dictPtr_;
  };

  void clear()
  {
    if(type_ == STRING || type_ == IDENTIFIER)
      {
      delete string_;
      }
  }

  void assignData(const vtkFoamToken& value)
  {
    switch(value.type_)
      {
      case PUNCTUATION:
        char_ = value.char_;
        break;
      case LABEL:
        int_ = value.int_;
        break;
      case SCALAR:
        double_ = value.double_;
        break;
      case STRING: case IDENTIFIER:
        string_ = new vtkStdString(*value.string_);
      }
  }

public:
  vtkFoamToken(): type_(UNDEFINED) {}
  vtkFoamToken(const vtkFoamToken& value): type_(value.type_)
    { assignData(value); }
  ~vtkFoamToken() { clear(); }

  const tokenType type() const { return type_; }

  template<typename T> bool is() const;
  template<typename T> T to() const;

  const vtkStdString toString() const { return *string_; }
  const vtkStdString toIdentifier() const { return *string_; }

  void setBad() { clear(); type_ = ERROR; }
  void setIdentifier(const vtkStdString& idString)
  { operator=(idString); type_ = IDENTIFIER; }

  void operator=(const char value)
  { clear(); type_ = PUNCTUATION; char_ = value; }
  void operator=(const int value)
  { clear(); type_ = LABEL; int_ = value; }
  void operator=(const double value)
  { clear(); type_ = SCALAR; double_ = value; }
  void operator=(const char *value)
  { clear(); type_ = STRING; string_ = new vtkStdString(value); }
  void operator=(const vtkStdString& value)
  { clear(); type_ = STRING; string_ = new vtkStdString(value); }
  void operator=(const vtkFoamToken& value)
  { clear(); type_ = value.type_; assignData(value); }
  bool operator==(const char value) const
  { return type_ == PUNCTUATION && char_ == value; }
  bool operator==(const int value) const
  { return type_ == LABEL && int_ == value; }
  bool operator==(const vtkStdString& value) const
  { return type_ == STRING && *string_ == value; }
  bool operator!=(const vtkStdString& value) const
  { return type_ != STRING || *string_ != value; }
  bool operator!=(const char value) const { return !operator==(value); }

  friend vtksys_ios::ostringstream& operator<<(vtksys_ios::ostringstream& str,
    const vtkFoamToken& value)
  {
    switch(value.type())
      {
      case ERROR:
        str << "badToken (an unexpected EOF?)";
        break;
      case PUNCTUATION:
        str << value.char_;
        break;
      case LABEL:
        str << value.int_;
        break;
      case SCALAR:
        str << value.double_;
        break;
      case STRING: case IDENTIFIER:
        str << *value.string_;
        break;
      }
    return str;
  }
};

template<> inline bool vtkOpenFOAMReaderPrivate::vtkFoamToken::is<int>() const
{ return type_ == LABEL; }
template<> inline bool vtkOpenFOAMReaderPrivate::vtkFoamToken::is<float>() const
{ return type_ == LABEL || type_ == SCALAR; }
template<> inline bool vtkOpenFOAMReaderPrivate::vtkFoamToken::is<double>()
  const
{ return type_ == SCALAR; }
template<> inline int vtkOpenFOAMReaderPrivate::vtkFoamToken::to() const
{ return int_; }
template<> inline float vtkOpenFOAMReaderPrivate::vtkFoamToken::to() const
{ return type_ == LABEL ? int_ : double_; }
template<> inline double vtkOpenFOAMReaderPrivate::vtkFoamToken::to() const
{ return type_ == LABEL ? int_ : double_; }

//-----------------------------------------------------------------------------
// class vtkFoamFileStack
// list of variables that have to be saved when a file is included.
struct vtkOpenFOAMReaderPrivate::vtkFoamFileStack
{
protected:
  vtkStdString fileName_;
  FILE *file_;
  bool isCompressed_;
  z_stream z_;
  int zStatus_;
  int lineNumber_;
#if VTK_FOAMFILE_RECOGNIZE_LINEHEAD
  bool wasNewline_;
#endif

  // buffer pointers. using raw pointers for performance reason.
  unsigned char *inbuf_;
  unsigned char *outbuf_;
  unsigned char *bufPtr_;
  unsigned char *bufEndPtr_;

  vtkFoamFileStack()
    : fileName_(), file_(NULL), isCompressed_(false), zStatus_(Z_OK),
    lineNumber_(0),
#if VTK_FOAMFILE_RECOGNIZE_LINEHEAD
    wasNewline_(true),
#endif
    inbuf_(NULL), outbuf_(NULL), bufPtr_(NULL), bufEndPtr_(NULL)
    { z_.zalloc = Z_NULL; z_.zfree = Z_NULL; z_.opaque = Z_NULL; }

  void reset()
  {
    // path_ = "";
    file_ = NULL;
    isCompressed_ = false;
    // zStatus_ = Z_OK;
    z_.zalloc = Z_NULL;
    z_.zfree = Z_NULL;
    z_.opaque = Z_NULL;
    // lineNumber_ = 0;
#if VTK_FOAMFILE_RECOGNIZE_LINEHEAD
    wasNewline_ = true;
#endif

    inbuf_ = NULL;
    outbuf_ = NULL;
    // bufPtr_ = NULL;
    // bufEndPtr_ = NULL;
  }

public:
  const vtkStdString& fileName() const { return fileName_; }
  int lineNumber() const { return lineNumber_; }
};

//-----------------------------------------------------------------------------
// class vtkFoamFile
// read and tokenize the input.
struct vtkOpenFOAMReaderPrivate::vtkFoamFile: public vtkFoamFileStack
{
public:
  // #inputMode values
  enum inputModes { INPUT_MODE_MERGE, INPUT_MODE_OVERWRITE, INPUT_MODE_ERROR };

private:
  bool is13Positions_;
  inputModes inputMode_;

  // inclusion handling
  vtkFoamFileStack *stack_[VTK_FOAMFILE_INCLUDE_STACK_SIZE];
  int stackI_;
  vtkStdString casePath_;

  // declare and define as private
  vtkFoamFile();
  bool inflateNext(unsigned char *buf, int requestSize);
  int nextTokenHead();
  // hacks to keep exception throwing / recursive codes out-of-line to make
  // putBack(), getc() and readExpecting() inline expandable
  void throwDuplicatedPutBackException();
  void throwUnexpectedEOFException();
  void throwUnexpectedNondigitCharExecption(const int c);
  void throwUnexpectedTokenException(const char, const int c);
  int readNext();

  void putBack(const int c)
  {
    if(--bufPtr_ < outbuf_)
      {
      throwDuplicatedPutBackException();
      }
    *bufPtr_ = c;
  }

  // get a character
  int getc() { return bufPtr_ == bufEndPtr_ ? readNext() : *bufPtr_++; }

  vtkFoamError stackString()
  {
    vtksys_ios::ostringstream os;
    if(stackI_ > 0)
      {
      os << "\n included";

      for(int stackI = stackI_ - 1; stackI >= 0; stackI--)
        {
        os << " from line " << stack_[stackI]->lineNumber() << " of "
          << stack_[stackI]->fileName() << "\n";
        }
      os << ": ";
      }
    return vtkFoamError() << os.str();
  }

  bool closeIncludedFile()
  {
    if(stackI_ == 0)
      {
      return false;
      }
    clear();
    stackI_--;
    // use the default bitwise assignment operator
    vtkFoamFileStack::operator=(*stack_[stackI_]);
    delete stack_[stackI_];
    return true;
  }

  void clear()
  {
    if(isCompressed_)
      {
      inflateEnd(&z_);
      }

    delete [] inbuf_;
    delete [] outbuf_;
    inbuf_ = outbuf_ = NULL;

    if(file_ != NULL)
      {
      fclose(file_);
      file_ = NULL;
      }
    // don't reset the line number so that the last line number is
    // retained after close
    // lineNumber_ = 0;
  }

  const vtkStdString extractPath(const vtkStdString& path) const
  {
#if defined(_WIN32)
    const vtkStdString pathFindSeparator = "/\\", pathSeparator = "\\";
#else
    const vtkStdString pathFindSeparator = "/", pathSeparator = "/";
#endif
    const vtkStdString::size_type pos = path.find_last_of(pathFindSeparator);
    return pos == vtkStdString::npos ? vtkStdString(".") + pathSeparator
      : path.substr(0, pos + 1);
  }

public:
  vtkFoamFile(const vtkStdString& casePath)
    : vtkFoamFileStack(), is13Positions_(false), inputMode_(INPUT_MODE_ERROR),
      stackI_(0), casePath_(casePath)
    {}
  ~vtkFoamFile() { close(); }

  void setIs13Positions(const bool is13Positions)
    { is13Positions_ = is13Positions; }
  const bool getIs13Positions() const { return is13Positions_; }
  inputModes getInputMode() const { return inputMode_; }
  const vtkStdString casePath() const { return casePath_; }
  const vtkStdString filePath() const { return extractPath(fileName_); }

  vtkStdString expandPath(const vtkStdString& pathIn,
    const vtkStdString& defaultPath)
  {
    vtkStdString expandedPath;
    bool isExpanded = false, wasPathSeparator = true;
    const int nChars = pathIn.length();
    for(int charI = 0; charI < nChars;)
      {
      char c = pathIn[charI];
      switch(c)
        {
        case '$': // $-variable expansion
          {
          vtkStdString variable;
          while(++charI < nChars && (isalnum(pathIn[charI])
            || pathIn[charI] == '_'))
            {
            variable += pathIn[charI];
            }
          if(variable == "FOAM_CASE") // discard path until the variable
            {
            expandedPath = casePath_;
	    wasPathSeparator = true;
	    isExpanded = true;
            }
          else
            {
            const char *value = getenv(variable.c_str());
            if(value != NULL)
              {
              expandedPath += value;
              }
            const vtkStdString::size_type len = expandedPath.length();
            if(len > 0)
              {
              const char c = expandedPath[len - 1];
              wasPathSeparator = (c == '/' || c == '\\');
              }
            else
              {
              wasPathSeparator = false;
              }
            }
          }
          break;
        case '~': // home directory expansion
	  // not using kwsys::SystemTools::ConvertToUnixSlashes() for
	  // a bit better handling of "~"
	  if(wasPathSeparator)
            {
            vtkStdString userName;
            while(++charI < nChars && (pathIn[charI] != '/'
              && pathIn[charI] != '\\') && pathIn[charI] != '$')
              {
              userName += pathIn[charI];
              }
            if(userName == "")
              {
	      const char *homePtr = getenv("HOME");
	      if(homePtr == NULL)
		{
#if defined(_WIN32) && !defined(__CYGWIN__) || defined(__LIBCATAMOUNT__)
		expandedPath = "";
#else
		const struct passwd *pwentry = getpwuid(getuid());
	        if(pwentry == NULL)
		  {
                  throw stackString() << "Home directory path not found";
		  }
                expandedPath = pwentry->pw_dir;
#endif
		}
	      else
		{
	        expandedPath = homePtr;
		}
              }
            else
              {
#if defined(_WIN32) && !defined(__CYGWIN__) || defined(__LIBCATAMOUNT__)
	      const char *homePtr = getenv("HOME");
              expandedPath = extractPath(homePtr ? homePtr : "") + userName;
#else
	      if(userName == "OpenFOAM")
		{
		// so far only "~/.OpenFOAM" expansion is supported
		const char *homePtr = getenv("HOME");
		if(homePtr == NULL)
		  {
		  expandedPath = "";
		  }
		else
		  {
		  expandedPath = vtkStdString(homePtr) + "/.OpenFOAM";
		  }
		}
	      else
		{
		const struct passwd *pwentry = getpwnam(userName.c_str());
	        if(pwentry == NULL)
		  {
		  throw stackString() << "Home directory for user "
                    << userName.c_str() << " not found";
		  }
                expandedPath = pwentry->pw_dir;
		}
#endif
              }
	    wasPathSeparator = false;
            isExpanded = true;
	    break;
            }
	  // fall through
        default:
	  wasPathSeparator = (c == '/' || c == '\\');
          expandedPath += c;
          charI++;
        }
      }
    if(isExpanded || expandedPath.substr(0, 1) == "/"
      || expandedPath.substr(0, 1) == "\\")
      {
      return expandedPath;
      }
    else
      {
      return defaultPath + expandedPath;
      }
  }

  void includeFile(const vtkStdString& includedFileName,
    const vtkStdString& defaultPath)
  {
    if(stackI_ >= VTK_FOAMFILE_INCLUDE_STACK_SIZE)
      {
      throw stackString() << "Exceeded maximum #include recursions of "
        << VTK_FOAMFILE_INCLUDE_STACK_SIZE;
      }
    // use the default bitwise copy constructor
    stack_[stackI_++] = new vtkFoamFileStack(*this);
    vtkFoamFileStack::reset();

    open(expandPath(includedFileName, defaultPath));
  }

  // the tokenizer
  // returns true if success, false if encountered EOF
  bool read(vtkFoamToken& token)
  {
    // expanded the outermost loop in nextTokenHead() for performance
    int c;
    while(isspace(c = getc())) // isspace() accepts -1 as EOF
      {
      if(c == '\n')
        {
        ++lineNumber_;
#if VTK_FOAMFILE_RECOGNIZE_LINEHEAD
        wasNewline_ = true;
#endif
        }
      }
    if(c == 47) // '/' == 47
      {
      putBack(c);
      c = nextTokenHead();
      }
#if VTK_FOAMFILE_RECOGNIZE_LINEHEAD
    if(c != '#')
      {
      wasNewline_ = false;
      }
#endif

    const int MAXLEN =  1024;
    char buf[MAXLEN + 1];
    int charI = 0;
    switch(c)
      {
      case '(': case ')':
        // high-priority punctuation token
        token = (char)c;
        return true;
      case '1': case '2': case '3': case '4': case '5': case '6': case '7':
      case '8': case '9': case '0': case '-':
        // undetermined number token
        do
          {
          buf[charI++] = c;
          } while(isdigit(c = getc()) && charI < MAXLEN);
        if(c != '.' && c != 'e' && c != 'E' && charI < MAXLEN && c != EOF)
          {
          // label token
          buf[charI] = '\0';
          token = static_cast<int>(strtol(buf, NULL, 10));
          putBack(c);
          return true;
          }
        // fall through
      case '.':
        // scalar token
        if(c == '.' && charI < MAXLEN)
          {
          // read decimal fraction part
          buf[charI++] = c;
          while(isdigit(c = getc()) && charI < MAXLEN)
            {
            buf[charI++] = c;
            }
          }
        if((c == 'e' || c == 'E') && charI < MAXLEN)
          {
          // read exponent part
          buf[charI++] = c;
          if(((c = getc()) == '+' || c == '-') && charI < MAXLEN)
            {
            buf[charI++] = c;
            c = getc();
            }
          while(isdigit(c) && charI < MAXLEN)
            {
            buf[charI++] = c;
            c = getc();
            }
          }
        if(charI == 1 && buf[0] == '-')
          {
          token = '-';
          putBack(c);
          return true;
          }
        buf[charI] = '\0';
        token = strtod(buf, NULL);
        putBack(c);
        break;
      case ';': case '{': case '}': case '[': case ']': case ':': case ',':
      case '=': case '+': case '*': case '/':
        // low-priority punctuation token
        token = (char)c;
        return true;
      case '"':
        {
        // string token
        bool wasEscape = false;
        while((c = getc()) != EOF && charI < MAXLEN)
          {
          if(c == '\\' && !wasEscape)
            {
            wasEscape = true;
            continue;
            }
          else if(c == '"' && !wasEscape)
            {
            break;
            }
          else if(c == '\n')
            {
            ++lineNumber_;
            if(!wasEscape)
              {
              throw stackString() << "Unescaped newline in string constant";
              }
            }
          buf[charI++] = c;
          wasEscape = false;
          }
        buf[charI] = '\0';
        token = buf;
        }
        break;
      case EOF:
        // end of file
        token.setBad();
        return false;
      case '$':
        {
        vtkFoamToken identifierToken;
        if(!read(identifierToken))
          {
          throw stackString() << "Unexpected EOF reading identifier";
          }
        if(identifierToken.type() != vtkFoamToken::STRING)
          {
          throw stackString() << "Expected a word, found " << identifierToken;
          }
        token.setIdentifier(identifierToken.toString());
        return true;
        }
      case '#':
        {
#if VTK_FOAMFILE_RECOGNIZE_LINEHEAD
        // placing #-directives in the middle of a line looks like
        // valid for the genuine OF 1.5 parser
        if(!wasNewline_)
          {
          throw stackString()
            << "Encountered #-directive in the middle of a line";
          }
        wasNewline_ = false;
#endif
        // read directive
        vtkFoamToken directiveToken;
        if(!read(directiveToken))
          {
          throw stackString() << "Unexpected EOF reading directive";
          }
        if(directiveToken == "include")
          {
          vtkFoamToken fileNameToken;
          if(!read(fileNameToken))
            {
            throw stackString() << "Unexpected EOF reading filename";
            }
          includeFile(fileNameToken.toString(), extractPath(fileName_));
          }
        else if(directiveToken == "inputMode")
          {
          vtkFoamToken modeToken;
          if(!read(modeToken))
            {
            throw stackString() << "Unexpected EOF reading inputMode specifier";
            }
          if(modeToken == "merge")
            {
            inputMode_ = INPUT_MODE_MERGE;
            }
          else if(modeToken == "overwrite")
            {
            inputMode_ = INPUT_MODE_OVERWRITE;
            }
          else if(modeToken == "error" || modeToken == "default")
            {
            inputMode_ = INPUT_MODE_ERROR;
            }
          else
            {
            throw stackString() << "Expected one of inputMode specifiers "
              "(merge, overwrite, error, default), found " << modeToken;
            }
          }
        else
          {
          throw stackString() << "Unsupported directive " << directiveToken;
          }
        return read(token);
        }
      default:
        // parses as a word token, but gives the STRING type for simplicity
        int inBrace = 0;
        do
          {
          if(c == '(')
            {
            inBrace++;
            }
          else if(c == ')' && --inBrace == -1)
            {
            break;
            }
          buf[charI++] = c;
          // valid characters that constitutes a word
          // cf. src/OpenFOAM/primitives/strings/word/wordI.H
          } while((c = getc()) != EOF && !isspace(c) && c != '"' && c != '/'
          && c != ';' && c != '{' && c != '}' && charI < MAXLEN);
        buf[charI] = '\0';
        token = buf;
        putBack(c);
      }

    if(c == EOF)
      {
      throwUnexpectedEOFException();
      }
    if(charI == MAXLEN)
      {
      throw stackString() << "Exceeded maximum allowed length of " << MAXLEN
        << " chars";
      }
    return true;
  }

  void open(const vtkStdString& fileName)
  {
    // reset line number to indicate the beginning of the file when an
    // exception is thrown
    lineNumber_ = 0;
    fileName_ = fileName;

    if(file_ != NULL)
      {
      throw stackString() << "File already opened within this object";
      }

    if((file_ = fopen(fileName.c_str(), "rb")) == NULL)
      {
      throw stackString() << "Can't open";
      }

    unsigned char zMagic[2];
    if(fread(zMagic, 1, 2, file_) == 2
      && zMagic[0] == 0x1f && zMagic[1] == 0x8b)
      {
      // gzip-compressed format
      z_.avail_in = 0;
      z_.next_in = Z_NULL;
      // + 32 to automatically recognize gzip format
      if(inflateInit2(&z_, 15 + 32) == Z_OK)
        {
        isCompressed_ = true;
        inbuf_ = new unsigned char[VTK_FOAMFILE_INBUFSIZE];
        }
      else
        {
        fclose(file_);
        file_ = NULL;
        throw stackString() << "Can't init zstream " << (z_.msg ? z_.msg : "");
        }
      }
    else
      {
      // uncompressed format
      isCompressed_ = false;
      }
    rewind(file_);

    zStatus_ = Z_OK;
    outbuf_ = new unsigned char[VTK_FOAMFILE_OUTBUFSIZE + 1];
    bufPtr_ = outbuf_ + 1;
    bufEndPtr_ = bufPtr_;
    lineNumber_ = 1;
  }

  void close()
  {
    while(closeIncludedFile());
    clear();
  }

  // gzread with buffering handling
  int read(unsigned char *buf, const int len)
  {
    int readlen;
    const int buflen = bufEndPtr_ - bufPtr_;
    if(len > buflen)
      {
      memcpy(buf, bufPtr_, buflen);
      readlen = inflateNext(buf + buflen, len - buflen);
      if(readlen >= 0)
        {
        readlen += buflen;
        }
      else
        {
        if(buflen == 0) // return EOF
          {
          readlen = -1;
          }
        else
          {
          readlen = buflen;
          }
        }
      bufPtr_ = bufEndPtr_;
      }
    else
      {
      memcpy(buf, bufPtr_, len);
      bufPtr_ += len;
      readlen = len;
      }
    for(int i = 0; i < readlen; i++)
      {
      if(buf[i] == '\n')
        {
        lineNumber_++;
        }
      }
    return readlen;
  }

  void readExpecting(const char expected)
  {
    // skip prepending invalid chars
    // expanded the outermost loop in nextTokenHead() for performance
    int c;
    while(isspace(c = getc())) // isspace() accepts -1 as EOF
      {
      if(c == '\n')
        {
        ++lineNumber_;
#if VTK_FOAMFILE_RECOGNIZE_LINEHEAD
        wasNewline_ = true;
#endif
        }
      }
    if(c == 47) // '/' == 47
      {
      putBack(c);
      c = nextTokenHead();
      }
    if(c != expected)
      {
      throwUnexpectedTokenException(expected, c);
      }
  }

  void readExpecting(const char* str)
  {
    vtkFoamToken t;
    if(!read(t) || t != str)
      {
      throw stackString() << "Expected string \"" << str << "\", found " << t;
      }
  }

  template<class T> T readValue(); //{ return static_cast<T>(0); }
};

int vtkOpenFOAMReaderPrivate::vtkFoamFile::readNext()
{
  if(!inflateNext(outbuf_ + 1, VTK_FOAMFILE_OUTBUFSIZE))
    {
    return closeIncludedFile() ? getc() : EOF;
    }
  return *bufPtr_++;
}

// specialized for reading an integer value.
// not using the standard strtol() for speed reason.
template<> int vtkOpenFOAMReaderPrivate::vtkFoamFile::readValue()
{
  // skip prepending invalid chars
  // expanded the outermost loop in nextTokenHead() for performance
  int c;
  while(isspace(c = getc())) // isspace() accepts -1 as EOF
    {
    if(c == '\n')
      {
      ++lineNumber_;
#if VTK_FOAMFILE_RECOGNIZE_LINEHEAD
      wasNewline_ = true;
#endif
      }
    }
  if(c == 47) // '/' == 47
    {
    putBack(c);
    c = nextTokenHead();
    }

  int nonNegative = c - 45; // '-' == 45
  if(nonNegative == 0 || c == 43) // '+' == 43
    {
    c = getc();
    if(c == '\n')
      {
      ++lineNumber_;
#if VTK_FOAMFILE_RECOGNIZE_LINEHEAD
      wasNewline_ = true;
#endif
      }
    }

  if(!isdigit(c)) // isdigit() accepts -1 as EOF
    {
    if(c == EOF)
      {
      throwUnexpectedEOFException();
      }
    else
      {
      throwUnexpectedNondigitCharExecption(c);
      }
    }

  int num = c - 48; // '0' == 48
  while(isdigit(c = getc()))
    {
    num = 10 * num + c - 48;
    }

  if(c == EOF)
    {
    throwUnexpectedEOFException();
    }
  putBack(c);

  return nonNegative ? num : -num;
}

// extreamely simplified high-performing string to floating point
// conversion code based on
// ParaView3/VTK/Utilities/vtksqlite/vtk_sqlite3.c
template<> float vtkOpenFOAMReaderPrivate::vtkFoamFile::readValue()
{
  // skip prepending invalid chars
  // expanded the outermost loop in nextTokenHead() for performance
  int c;
  while(isspace(c = getc())) // isspace() accepts -1 as EOF
    {
    if(c == '\n')
      {
      ++lineNumber_;
#if VTK_FOAMFILE_RECOGNIZE_LINEHEAD
      wasNewline_ = true;
#endif
      }
    }
  if(c == 47) // '/' == 47
    {
    putBack(c);
    c = nextTokenHead();
    }

  // determine sign
  int nonNegative = c - 45; // '-' == 45
  if(nonNegative == 0 || c == 43) // '+' == 43
    {
    c = getc();
    if(c == '\n')
      {
      ++lineNumber_;
#if VTK_FOAMFILE_RECOGNIZE_LINEHEAD
      wasNewline_ = true;
#endif
      }
    }

  if(!isdigit(c) && c != 46) // '.' == 46, isdigit() accepts EOF
    {
    throwUnexpectedNondigitCharExecption(c);
    }

  // read integer part
  double num = c - 48; // '0' == 48
  while(isdigit(c = getc()))
    {
    num = num * 10.0 + (c - 48);
    }

  // read decimal part
  if(c == 46) // '.'
    {
    double divisor = 1.0;

    while(isdigit(c = getc()))
      {
      num = num * 10.0 + (c - 48);
      divisor *= 10.0;
      }
    num /= divisor;
    }

  // read exponent part
  if(c == 69 || c == 101) // 'E' == 69, 'e' == 101
    {
    int esign = 1;
    int eval = 0;
    double scale = 1.0;

    c = getc();
    if(c == 45) // '-'
      {
      esign = -1;
      c = getc();
      }
    else if(c == 43) // '+'
      {
      c = getc();
      }

    while(isdigit(c))
      {
      eval = eval * 10 + (c - 48);
      c = getc();
      }

    // fast exponent multiplication!
    while(eval >= 64)
      {
      scale *= 1.0e+64;
      eval -= 64;
      }
    while(eval >= 16)
      {
      scale *= 1.0e+16;
      eval -= 16;
      }
    while(eval >= 4)
      {
      scale *= 1.0e+4;
      eval -= 4;
      }
    while(eval >= 1)
      {
      scale *= 1.0e+1;
      eval -= 1;
      }

    if(esign < 0)
      {
      num /= scale;
      }
    else
      {
      num *= scale;
      }
    }

  if(c == EOF)
    {
    throwUnexpectedEOFException();
    }
  putBack(c);

  return static_cast<float>(nonNegative ? num : -num);
}

// hacks to keep exception throwing code out-of-line to make
// putBack() and readExpecting() inline expandable
void vtkOpenFOAMReaderPrivate::vtkFoamFile::throwUnexpectedEOFException()
{ throw stackString() << "Unexpected EOF"; }

void vtkOpenFOAMReaderPrivate::vtkFoamFile
  ::throwUnexpectedNondigitCharExecption(const int c)
{
  throw stackString() << "Expected a number, found a non-digit character "
    << static_cast<char>(c);
}

void vtkOpenFOAMReaderPrivate::vtkFoamFile::throwUnexpectedTokenException(
  const char expected, const int c)
{
  vtkFoamError sstr;
  sstr << stackString() << "Expected punctuation token '" << expected
    << "', found ";
  if(c == EOF)
    {
    sstr << "EOF";
    }
  else
    {
    sstr << static_cast<char>(c);
    }
  throw sstr;
}

void vtkOpenFOAMReaderPrivate::vtkFoamFile::throwDuplicatedPutBackException()
{ throw stackString() << "Attempted duplicated putBack()"; }

bool vtkOpenFOAMReaderPrivate::vtkFoamFile::inflateNext(unsigned char *buf,
  int requestSize)
{
  int size;
  if(isCompressed_)
    {
    if(zStatus_ != Z_OK)
      {
      return false;
      }
    z_.next_out = buf;
    z_.avail_out = requestSize;

    do
      {
      if(z_.avail_in == 0)
        {
        z_.next_in = inbuf_;
        z_.avail_in = fread(inbuf_, 1, VTK_FOAMFILE_INBUFSIZE, file_);
        if(ferror(file_))
          {
          throw stackString() << "Fread failed";
          }
        }
      zStatus_ = inflate(&z_, Z_NO_FLUSH);
      if(zStatus_ == Z_STREAM_END
#if VTK_FOAMFILE_OMIT_CRCCHECK
        // the dummy CRC function causes data error when finalizing
        // so we have to proceed even when a data error is detected
        || zStatus_ == Z_DATA_ERROR
#endif
        )
        {
        break;
        }
      if(zStatus_ != Z_OK)
        {
        throw stackString() << "Inflation failed: " << (z_.msg ? z_.msg : "");
        }
      }
    while(z_.avail_out > 0);

    size = requestSize - z_.avail_out;
    }
  else
    {
    // not compressed
    size = fread(buf, 1, requestSize, file_);
    }

  if(size <= 0)
    {
    // retain the current location bufPtr_ to the end of the buffer so that
    // getc() returns EOF again when called next time
    return false;
    }
  // size > 0
  bufPtr_ = outbuf_ + 1; // reserve the first byte for getback char
  bufEndPtr_ = bufPtr_ + size;
  return true;
}

// get next semantically valid character
int vtkOpenFOAMReaderPrivate::vtkFoamFile::nextTokenHead()
{
  for(;;)
    {
    int c;
    while(isspace(c = getc())) // isspace() accepts -1 as EOF
      {
      if(c == '\n')
        {
        ++lineNumber_;
#if VTK_FOAMFILE_RECOGNIZE_LINEHEAD
        wasNewline_ = true;
#endif
        }
      }
    if(c == '/')
      {
      if((c = getc()) == '/')
        {
        while((c = getc()) != EOF && c != '\n');
        if(c == EOF)
          {
          return c;
          }
        ++lineNumber_;
#if VTK_FOAMFILE_RECOGNIZE_LINEHEAD
        wasNewline_ = true;
#endif
        }
      else if(c == '*')
        {
        for(;;)
          {
          while((c = getc()) != EOF && c != '*')
            {
            if(c == '\n')
              {
              ++lineNumber_;
              }
            }
          if(c == EOF)
            {
            return c;
            }
          else if((c = getc()) == '/')
            {
            break;
            }
          putBack(c);
          }
        }
      else
        {
        putBack(c); // may be an EOF
        return '/';
        }
      }
    else // may be an EOF
      {
      return c;
      }
    }
}

//-----------------------------------------------------------------------------
// class vtkFoamIOobject
// holds file handle, file format, name of the object the file holds and
// type of the object.
struct vtkOpenFOAMReaderPrivate::vtkFoamIOobject: public vtkFoamFile
{
public:
  enum fileFormat { UNDEFINED, ASCII, BINARY };

private:
  fileFormat format_;
  vtkStdString objectName_;
  vtkStdString headerClassName_;
  vtkFoamError e_;

  vtkFoamIOobject();
  void readHeader(); // defined later
public:
  vtkFoamIOobject(const vtkStdString& casePath)
    : vtkFoamFile(casePath), format_(UNDEFINED), e_() {}
  ~vtkFoamIOobject() { close(); }

  bool open(const vtkStdString& file)
  {
    try
      {
      vtkFoamFile::open(file);
      }
    catch(vtkFoamError& e)
      {
      e_ = e;
      return false;
      }

    try
      {
      readHeader();
      }
    catch(vtkFoamError& e)
      {
      vtkFoamFile::close();
      e_ = e;
      return false;
      }
    return true;
  }

  void close()
  {
    vtkFoamFile::close();
    format_ = UNDEFINED;
    objectName_.erase();
    headerClassName_.erase();
    e_.erase();
  }
  const fileFormat format() const { return format_; }
  const vtkStdString& className() const { return headerClassName_; }
  const vtkStdString& objectName() const { return objectName_; }
  const vtkFoamError& error() const { return e_; }
  void setError(const vtkFoamError& e) { e_ = e; }
};

//-----------------------------------------------------------------------------
// class vtkFoamEntryValue
// a class that represents a value of a dictionary entry that corresponds to
// its keyword. note that an entry can have more than one value.
struct vtkOpenFOAMReaderPrivate::vtkFoamEntryValue: public vtkFoamToken
{
private:
  bool isUniform_;
  bool managed_;
  vtkFoamEntry *upperEntryPtr_;

  vtkFoamEntryValue();
  vtkObjectBase *toVTKObject() { return vtkObjectPtr_; }
  void clear();
  void readList(vtkFoamIOobject& io);

public:
  // reads primitive int/float lists
  template <typename listT, typename primitiveT> class listTraits
  {
    listT *ptr_;

  public:
    listTraits(): ptr_(listT::New()) {}
    listT *ptr() { return ptr_; }
    void readUniformValues(vtkFoamIOobject& io, const int size)
    {
      primitiveT value = io.readValue<primitiveT>();
      for(int i = 0; i < size; i++)
	{
	ptr_->SetValue(i, value);
	}
    }
    void readAsciiList(vtkFoamIOobject& io, const int size)
    {
      for(int i = 0; i < size; i++)
        {
	ptr_->SetValue(i, io.readValue<primitiveT>());
        }
    }
    void readBinaryList(vtkFoamIOobject& io, const int size)
    {
      io.read(reinterpret_cast<unsigned char *>(ptr_->GetPointer(0)),
        size * sizeof(primitiveT));
    }
    void readValue(vtkFoamIOobject&, vtkFoamToken& currToken)
    {
      if(!currToken.is<primitiveT>())
        {
        throw vtkFoamError() << "Expected an integer or a (, found "
          << currToken;
        }
      ptr_->InsertNextValue(currToken.to<primitiveT>());
    }
  };

  // reads rank 1 lists of types vector, sphericalTensor, symmTensor
  // and tensor. if isPositions is true it reads Cloud type of data as
  // particle positions. cf. (the positions format)
  // src/lagrangian/basic/particle/particleIO.C
  template <typename listT, typename primitiveT,
    int nComponents, bool isPositions> class vectorListTraits
  {
    listT *ptr_;

  public:
    vectorListTraits(): ptr_(listT::New())
    { ptr_->SetNumberOfComponents(nComponents); }
    listT *ptr() { return ptr_; }
    void readUniformValues(vtkFoamIOobject& io, const int size)
    {
      io.readExpecting('(');
      primitiveT vectorValue[nComponents];
      for(int j = 0; j < nComponents; j++)
        {
	vectorValue[j] = io.readValue<primitiveT>();
        }
      for(int i = 0; i < size; i++)
        {
        ptr_->SetTuple(i, vectorValue);
        }
      io.readExpecting(')');
      if(isPositions)
        {
        // skip label celli
	io.readValue<int>();
        }
    }
    void readAsciiList(vtkFoamIOobject& io, const int size)
    {
      for(int i = 0; i < size; i++)
        {
        io.readExpecting('(');
        primitiveT *vectorTupleI = ptr_->GetPointer(nComponents * i);
        for(int j = 0; j < nComponents; j++)
          {
          vectorTupleI[j] = io.readValue<primitiveT>();
          }
        io.readExpecting(')');
        if(isPositions)
          {
          // skip label celli
	  io.readValue<int>();
          }
        }
    }
    void readBinaryList(vtkFoamIOobject& io, const int size)
    {
      if(isPositions) // lagrangian/positions (class Cloud)
        {
        // allocate space along with the larger 1.4 format since the
        // size must be determined at compile-time. we allocate on the
        // stack to avoid leak when an exception is thrown.
        unsigned char buffer[sizeof(double) * (nComponents + 1)
	  + 2 * sizeof(int)];
        const int nBytes = (io.getIs13Positions()
        // skip label celli
          ? sizeof(double) * nComponents + sizeof(int)
        // skip label celli, label facei and scalar stepFraction
          : sizeof(double) * (nComponents + 1) + 2 * sizeof(int));
        for(int i = 0; i < size; i++)
          {
          io.readExpecting('(');
          io.read(buffer, nBytes);
          ptr_->SetTuple(i, reinterpret_cast<double *>(buffer));
          io.readExpecting(')');
          }
	}
      else
	{
        for(int i = 0; i < size; i++)
          {
          double buffer[nComponents];
          io.read(reinterpret_cast<unsigned char *>(buffer),
	    sizeof(double) * nComponents);
          ptr_->SetTuple(i, buffer);
          }
	}
    }
    void readValue(vtkFoamIOobject& io, vtkFoamToken& currToken)
    {
      if(currToken != '(')
        {
        throw vtkFoamError() << "Expected '(', found " << currToken;
        }
      primitiveT v[nComponents];
      for(int j = 0; j < nComponents; j++)
        {
	v[j] = io.readValue<primitiveT>();
        }
      ptr_->InsertNextTuple(v);
      io.readExpecting(')');
    }
  };

  vtkFoamEntryValue(vtkFoamEntry *upperEntryPtr)
    : vtkFoamToken(), isUniform_(false), managed_(true),
      upperEntryPtr_(upperEntryPtr) {}
  vtkFoamEntryValue(vtkFoamEntryValue&, vtkFoamEntry *);
  ~vtkFoamEntryValue() { clear(); }

  void setEmptyList() { clear(); isUniform_ = false; type_ = EMPTYLIST; }
  bool isUniform() const { return isUniform_; }
  void read(vtkFoamIOobject& io);
  void readDictionary(vtkFoamIOobject& io, const vtkFoamToken& firstKeyword);
  vtkIntArray& labelList() const { return *labelListPtr_; }
  const intVectorVector& labelListList() const { return *labelListListPtr_; }
  vtkFloatArray& scalarList() { return *scalarListPtr_; }
  vtkFloatArray& vectorList() { return *vectorListPtr_; }
  vtkFoamDict& dictionary() { return *dictPtr_; }

  void *ptr()
  {
    managed_ = false; // the returned pointer will not be deleted by the d'tor
    return (void *)labelListPtr_; // all list pointers are in a single union
  }

  vtkStdString toString() const
  { return type_ == STRING ? vtkFoamToken::toString() : vtkStdString(); }
  float toFloat() const
  {
    return type_ == SCALAR || type_ == LABEL ? vtkFoamToken::to<float>() : 0.0F;
  }
  double toDouble() const
  {
    return type_ == SCALAR || type_ == LABEL ? vtkFoamToken::to<double>() : 0.0;
  }
  int toInt() const
  { return type_ == LABEL ? vtkFoamToken::to<int>() : 0; }

  // the following two are for an exceptional expression of
  // `LABEL{LABELorSCALAR}' without type prefix (e. g. `2{-0}' in
  // mixedRhoE B.C. in rhopSonicFoam/shockTube)
  void makeLabelList(const int labelValue, const int size)
  {
    labelListPtr_ = vtkIntArray::New();
    type_ = LABELLIST;
    vtkIntArray &ll = *labelListPtr_;
    ll.SetNumberOfValues(size);
    for(int i = 0; i < size; i++)
      {
      ll.SetValue(i, labelValue);
      }
  }
  void makeScalarList(const float scalarValue, const int size)
  {
    scalarListPtr_ = vtkFloatArray::New();
    type_ = SCALARLIST;
    vtkFloatArray& sl = *scalarListPtr_;
    sl.SetNumberOfValues(size);
    for(int i = 0; i < size; i++)
      {
      sl.SetValue(i, scalarValue);
      }
  }

  // reads dimensionSet
  void readDimensionSet(vtkFoamIOobject& io)
  {
    const int nDims = 7;
    labelListPtr_ = vtkIntArray::New();
    type_ = LABELLIST;
    labelListPtr_->SetNumberOfValues(nDims);
    for(int dimI = 0; dimI < nDims; dimI++)
      {
      labelListPtr_->SetValue(dimI, io.readValue<int>());
      }
    io.readExpecting(']');
  }

  // generic reader for nonuniform lists. requires size prefix of the
  // list to be present in the stream if the format is binary.
  template <vtkFoamToken::tokenType listType, typename traitsT>
  void readNonuniformList(vtkFoamIOobject& io)
  {
    vtkFoamToken currToken;
    if(!io.read(currToken))
      {
      throw vtkFoamError() << "Unexpected EOF";
      }
    traitsT list;
    type_ = listType;
    vtkObjectPtr_ = list.ptr();
    if(currToken.is<int>())
      {
      const int size = currToken.to<int>();
      if(size < 0)
        {
        throw vtkFoamError() << "List size must not be negative: size = "
          << size;
        }
      list.ptr()->SetNumberOfTuples(size);
      if(!io.read(currToken))
        {
        throw vtkFoamError() << "Unexpected EOF";
        }
      // some objects have lists with only one element enclosed by {}
      // e. g. simpleFoam/pitzDaily3Blocks/constant/polyMesh/faceZones
      if(currToken == '{')
        {
	list.readUniformValues(io, size);
        io.readExpecting('}');
	return;
        }
      else if(currToken != '(')
        {
        throw vtkFoamError() << "Expected '(', found " << currToken;
        }
      else if(io.format() == vtkFoamIOobject::ASCII)
        {
	list.readAsciiList(io, size);
        }
      else
        {
        if(size > 0) // avoid invalid access to ll.at(0)
          {
	  list.readBinaryList(io, size);
          }
        }
      io.readExpecting(')');
      }
    else if(currToken == '(')
      {
      while(io.read(currToken) && currToken != ')')
        {
	list.readValue(io, currToken);
        }
      list.ptr()->Squeeze();
      }
    else
      {
      throw vtkFoamError() << "Expected integer or '(', found " << currToken;
      }
  }

  // reads a list of labelLists. requires size prefix of the listList
  // to be present. size of each sublist must also be present in the
  // stream if the format is binary.
  void readLabelListList(vtkFoamIOobject& io)
  {
    vtkFoamToken currToken;
    if(!io.read(currToken))
      {
      throw vtkFoamError() << "Unexpected EOF";
      }
    if(currToken.type() == vtkFoamToken::LABEL)
      {
      const int sizeI = currToken.to<int>();
      if(sizeI < 0)
        {
        throw vtkFoamError() << "List size must not be negative: size = "
          << sizeI;
        }
      // gives initial guess for list size
      labelListListPtr_ = new intVectorVector(sizeI, 4 * sizeI);
      type_ = LABELLISTLIST;
      io.readExpecting('(');
      int bodyI = 0;
      for(int i = 0; i < sizeI; i++)
        {
        if(!io.read(currToken))
          {
          throw vtkFoamError() << "Unexpected EOF";
          }
        if(currToken.type() == vtkFoamToken::LABEL)
          {
	  const int sizeJ = currToken.to<int>();
          if(sizeJ < 0)
            {
            throw vtkFoamError() << "List size must not be negative: size = "
              << sizeJ;
            }
          if(bodyI + sizeJ > labelListListPtr_->bodySize())
            {
            const int newSize = labelListListPtr_->bodySize() + sizeJ;
            labelListListPtr_->resizeBody(newSize);
            }
          io.readExpecting('(');
          int *listI = labelListListPtr_->setIndex(i, bodyI);
          if(io.format() == vtkFoamIOobject::ASCII)
            {
            for(int j = 0; j < sizeJ; j++)
              {
	      listI[j] = io.readValue<int>();
              }
            }
          else
            {
            if(sizeJ > 0) // avoid invalid reference to labelListI.at(0)
              {
              io.read(reinterpret_cast<unsigned char*>(listI),
                sizeJ * sizeof(int));
              }
            }
          bodyI += sizeJ;
          io.readExpecting(')');
          }
        else if(currToken == '(')
          {
          labelListListPtr_->setIndex(i, bodyI);
          while(io.read(currToken) && currToken != ')')
            {
            if(currToken.type() != vtkFoamToken::LABEL)
              {
              throw vtkFoamError() << "Expected an integer, found "
                << currToken;
              }
            if(bodyI >= labelListListPtr_->bodySize())
              {
              const int newSize = labelListListPtr_->bodySize() + 1;
              labelListListPtr_->resizeBody(newSize);
              }
            labelListListPtr_->setValue(bodyI++, currToken.to<int>());
            }
          }
        else
          {
          throw vtkFoamError() << "Expected integer or '(', found "
            << currToken;
          }
        }
      // set the next index of the last element to calculate the last
      // subarray size
      labelListListPtr_->setIndex(sizeI, bodyI);
      // shrink to the actually used size
      labelListListPtr_->resizeBody(bodyI);
      io.readExpecting(')');
      }
    else
      {
      throw vtkFoamError() << "Expected integer, found " << currToken;
      }
  }

  bool readField(vtkFoamIOobject& io)
  {
    try
      {
      if(io.className() == "scalarField") // lagrangian scalars
        {
        readNonuniformList<SCALARLIST, listTraits<vtkFloatArray, float> >(
          io);
        }
      // polyMesh/points, lagrangian vectors
      else if(io.className() == "sphericalTensorField")
        {
        readNonuniformList<VECTORLIST,
          vectorListTraits<vtkFloatArray, float, 1, false> >(io);
        }
      else if(io.className() == "vectorField")
        {
        readNonuniformList<VECTORLIST,
          vectorListTraits<vtkFloatArray, float, 3, false> >(io);
        }
      else if(io.className() == "symmTensorField")
        {
        readNonuniformList<VECTORLIST,
          vectorListTraits<vtkFloatArray, float, 6, false> >(io);
        }
      else if(io.className() == "tensorField")
        {
        readNonuniformList<VECTORLIST,
          vectorListTraits<vtkFloatArray, float, 9, false> >(io);
        }
      else
        {
        throw vtkFoamError() << "Non-supported field type " << io.className();
        }
      }
    catch(vtkFoamError& e)
      {
      io.setError(e);
      return false;
      }
    return true;
  }
};

template<>
void vtkOpenFOAMReaderPrivate::vtkFoamEntryValue::listTraits<vtkFloatArray,
  float>::readBinaryList(vtkFoamIOobject& io, const int size)
{
  for(int i = 0; i < size; i++)
    {
    double buffer;
    io.read(reinterpret_cast<unsigned char *>(&buffer), sizeof(double));
    ptr_->SetValue(i, static_cast<float>(buffer));
    }
}

//-----------------------------------------------------------------------------
// class vtkFoamEntry
// a class that represents an entry of a dictionary. note that an
// entry can have more than one value.
struct vtkOpenFOAMReaderPrivate::vtkFoamEntry
  : public vtkstd::vector<vtkFoamEntryValue*>
{
private:
  typedef vtkstd::vector<vtkFoamEntryValue*> Superclass;
  vtkStdString keyword_;
  vtkFoamDict *upperDictPtr_;

  vtkFoamEntry();

public:
  vtkFoamEntry(vtkFoamDict *upperDictPtr): upperDictPtr_(upperDictPtr) {}
  vtkFoamEntry(const vtkFoamEntry& entry, vtkFoamDict *upperDictPtr)
    : Superclass(entry.size()), keyword_(entry.keyword()),
      upperDictPtr_(upperDictPtr)
  {
    for(size_t valueI = 0; valueI < entry.size(); valueI++)
      {
      Superclass::operator[](valueI)
	= new vtkFoamEntryValue(entry.value(valueI), this);
      }
  }

  ~vtkFoamEntry() { clear(); }

  void clear()
  {
    for(size_t i = 0; i < Superclass::size(); i++)
      {
      delete Superclass::operator[](i);
      }
    Superclass::clear();
  }
  const vtkStdString& keyword() const { return keyword_; }
  vtkStdString& keyword() { return keyword_; }
  size_t size() const { return Superclass::size(); }
  // returns false if the number of the values is 0 to simplify things
  bool found() const { return Superclass::size() > 0; }
  vtkFoamEntryValue& value(const int i) const
  { return *Superclass::operator[](i); }
  vtkFoamEntryValue& firstValue() const { return *Superclass::operator[](0); }
  vtkIntArray& labelList() const { return firstValue().labelList(); }
  const intVectorVector& labelListList() const
  { return firstValue().labelListList(); }
  vtkFloatArray& scalarList() { return firstValue().scalarList(); }
  vtkFloatArray& vectorList() { return firstValue().vectorList(); }
  vtkFoamDict& dictionary() // not using firstValue() for breaking constness
  { return Superclass::operator[](0)->dictionary(); }
  void *ptr() { return firstValue().ptr(); }
  vtkFoamDict *upperDictPtr() { return upperDictPtr_; }

  vtkStdString toString() const
  { return found() ? firstValue().toString() : vtkStdString(); }
  float toFloat() const { return found() ? firstValue().toFloat() : 0.0F; }
  double toDouble() const { return found() ? firstValue().toDouble() : 0.0; }
  int toInt() const { return found() ? firstValue().toInt() : 0; }

  void readDictionary(vtkFoamIOobject& io)
  {
    Superclass::push_back(new vtkFoamEntryValue(this));
    Superclass::back()->readDictionary(io, vtkFoamToken());
  }

  // read values of an entry
  void read(vtkFoamIOobject& io);
};

//-----------------------------------------------------------------------------
// class vtkFoamDict
// a class that holds a FoamFile data structure
struct vtkOpenFOAMReaderPrivate::vtkFoamDict: public vtkFoamEntryValue
{
private:
  vtkstd::vector<vtkFoamEntry*> entryPtrs_;
  vtkFoamEntry *dummyEntryPtr_;
  vtkFoamDict *upperDictPtr_;

  vtkFoamDict(const vtkFoamDict &);

public:
  vtkFoamDict(vtkFoamDict *upperDictPtr = NULL)
    : vtkFoamEntryValue(NULL), dummyEntryPtr_(NULL), upperDictPtr_(upperDictPtr)
  { dictPtr_  = NULL; } // avoid destruction in vtkFoamEntryValue::clear()
  vtkFoamDict(vtkFoamDict& dict, vtkFoamDict *upperDictPtr)
    : vtkFoamEntryValue(dict, NULL), entryPtrs_(dict.size()),
      dummyEntryPtr_(NULL), upperDictPtr_(upperDictPtr)
  {
    if(dict.type() == vtkFoamToken::DICTIONARY)
      {
      for(size_t entryI = 0; entryI < dict.size(); entryI++)
        {
        entryPtrs_[entryI] = new vtkFoamEntry(dict.entry(entryI), this);
        }
      }
  }

  ~vtkFoamDict()
  {
    if(type_ == DICTIONARY)
      {
      for(size_t i = 0; i < entryPtrs_.size(); i++)
        {
        delete entryPtrs_[i];
        }
      }
    delete dummyEntryPtr_;
  }

  size_t size() const { return entryPtrs_.size(); }
  vtkFoamDict *upperDictPtr() { return upperDictPtr_; }
  const vtkFoamEntry& entry(const int i) const { return *entryPtrs_[i]; }
  vtkFoamEntry& entry(const int i) { return *entryPtrs_[i]; }
  vtkFoamEntry& lookup(const vtkStdString& keyword)
  {
    if(type_ == DICTIONARY)
      {
      for(size_t i = 0; i < entryPtrs_.size(); i++)
        {
        if(entryPtrs_[i]->keyword() == keyword) // found
          {
          return *entryPtrs_[i];
          }
        }
      }

    // not found
    if(dummyEntryPtr_ == NULL)
      {
      dummyEntryPtr_ = new vtkFoamEntry(NULL);
      }
    return *dummyEntryPtr_;
  }

  // reads a FoamFile or a subdictionary. if the stream to be read is
  // a subdictionary the preceding '{' is assumed to have already been
  // thrown away.
  bool read(vtkFoamIOobject& io, const bool isSubDictionary = false,
    const vtkFoamToken& firstToken = vtkFoamToken())
  {
    try
      {
      vtkFoamToken currToken;
      if(firstToken.type() == vtkFoamToken::UNDEFINED)
        {
        // read the first token
        if(!io.read(currToken))
          {
          throw vtkFoamError() << "Unexpected EOF";
          }

        if(isSubDictionary)
          {
          // the following two else-if clauses are for an exceptional
          // expression of `LABEL{LABELorSCALAR}' without type prefix
          // (e. g. `2{-0}' in mixedRhoE B.C. in
          // rhopSonicFoam/shockTube)
          if(currToken.type() == vtkFoamToken::LABEL)
            {
            vtkFoamToken::operator=(currToken.to<int>());
            type_ = UNIFORMLABELLIST;
            io.readExpecting('}');
            return true;
            }
          else if(currToken.type() == vtkFoamToken::SCALAR)
            {
            vtkFoamToken::operator=(currToken.to<float>());
            type_ = UNIFORMSCALARLIST;
            io.readExpecting('}');
            return true;
            }
          // return as empty dictionary
          else if(currToken == '}')
            {
            type_ = DICTIONARY;
            return true;
            }
          }
        else
          {
          // list of dictionaries is read as a usual dictionary
          // polyMesh/boundary, point/face/cell-Zones
          if(currToken.type() == vtkFoamToken::LABEL)
            {
            io.readExpecting('(');
            if(currToken.to<int>() > 0)
              {
              if(!io.read(currToken))
                {
                throw vtkFoamError() << "Unexpected EOF";
                }
              // continue to read as a usual dictionary
              }
            else // return as empty dictionary
              {
              io.readExpecting(')');
              type_ = DICTIONARY;
              return true;
              }
            }
          // some boundary files does not have the number of boundary
          // patches (e.g. settlingFoam/tank3D). in this case we need to
          // explicitly read the file as a dictionary.
          else if(currToken == '('
            && io.className() == "polyBoundaryMesh") // polyMesh/boundary
            {
	    if(!io.read(currToken)) // read the first keyword
              {
              throw vtkFoamError() << "Unexpected EOF";
              }
            if(currToken == ')') // return as empty dictionary
              {
              type_ = DICTIONARY;
              return true;
              }
            }
          }
        }
      // if firstToken is given as string read the following stream as
      // subdictionary
      else if(firstToken.type() == vtkFoamToken::STRING)
        {
        type_ = DICTIONARY;
        entryPtrs_.push_back(new vtkFoamEntry(this));
        entryPtrs_.back()->keyword() = firstToken.toString();
        entryPtrs_.back()->readDictionary(io);
        if(!io.read(currToken) || currToken == '}' || currToken == ')')
          {
          return true;
          }
        }
      else // quite likely an identifier
        {
        currToken = firstToken;
        }

      if(currToken == ';' || currToken.type() == vtkFoamToken::STRING
        || currToken.type() == vtkFoamToken::IDENTIFIER) // general dictionary
        {
        // type must be set first so that lookup() works
        type_ = DICTIONARY;
        do
          {
          if(currToken.type() == vtkFoamToken::STRING)
            {
            vtkFoamEntry& previousEntry = lookup(currToken.toString());
            if(previousEntry.found())
              {
              if(io.getInputMode() == vtkFoamFile::INPUT_MODE_MERGE)
                {
                if(previousEntry.firstValue().type()
                  == vtkFoamToken::DICTIONARY)
                  {
                  io.readExpecting('{');
                  previousEntry.firstValue().dictionary().read(io, true);
                  }
                else
                  {
                  previousEntry.clear();
                  previousEntry.read(io);
                  }
                }
              else if(io.getInputMode() == vtkFoamFile::INPUT_MODE_OVERWRITE)
                {
                previousEntry.clear();
                previousEntry.read(io);
                }
              else // INPUT_MODE_ERROR
                {
                throw vtkFoamError() << "Found duplicated entries with keyword "
                  << currToken.toString();
                }
              }
            else
              {
              entryPtrs_.push_back(new vtkFoamEntry(this));
              entryPtrs_.back()->keyword() = currToken.toString();
              entryPtrs_.back()->read(io);
              }

            if(currToken == "FoamFile")
              {
              // delete the FoamFile header subdictionary entry
              delete entryPtrs_.back();
              entryPtrs_.pop_back();
              }
	    else if(currToken == "include")
              {
              // include the named file. Exiting the included file at
              // EOF will be handled automatically by
              // vtkFoamFile::closeIncludedFile()
              if(entryPtrs_.back()->firstValue().type() != vtkFoamToken::STRING)
                {
                throw vtkFoamError()
                  << "Expected string as the file name to be included, found "
                  << entryPtrs_.back()->firstValue();
                }
              vtkStdString includeFileName(entryPtrs_.back()->toString());
              delete entryPtrs_.back();
              entryPtrs_.pop_back();
              io.includeFile(includeFileName, io.filePath());
              }
            }
          else if(currToken.type() == vtkFoamToken::IDENTIFIER)
            {
	    // substitute identifier
	    const vtkStdString identifier(currToken.toIdentifier());

	    for(vtkFoamDict *uDictPtr = this;;)
	      {
	      const vtkFoamEntry&
		identifiedEntry = uDictPtr->lookup(identifier);

	      if(identifiedEntry.found())
		{
		if(identifiedEntry.firstValue().type()
                  != vtkFoamToken::DICTIONARY)
		  {
		  throw vtkFoamError()
		    << "Expected dictionary for substituting entry "
		    << identifier;
		  }
		vtkFoamDict& identifiedDict
		  = identifiedEntry.firstValue().dictionary();
		for(size_t entryI = 0; entryI < identifiedDict.size(); entryI++)
		  {
                  // I think #inputMode handling should be done here
                  // as well, but the genuine FoamFile parser for OF
                  // 1.5 does not seem to be doing it.
		  entryPtrs_.push_back(
                    new vtkFoamEntry(identifiedDict.entry(entryI), this));
		  }
                break;
                }
              else
                {
                uDictPtr = uDictPtr->upperDictPtr();
                if(uDictPtr == NULL)
                  {
                  throw vtkFoamError() << "Substituting entry " << identifier
                    << " not found";
                  }
                }
              }
            }
          // skip empty entry only with ';'
          } while(io.read(currToken)
          && (currToken.type() == vtkFoamToken::STRING
          || currToken.type() == vtkFoamToken::IDENTIFIER || currToken == ';'));

        if(currToken.type() == vtkFoamToken::ERROR || currToken == '}'
          || currToken == ')')
          {
          return true;
          }
	type_ = UNDEFINED;
        throw vtkFoamError()
          << "Expected keyword, closing brace, ';' or EOF, found " << currToken;
        }
      throw vtkFoamError() << "Expected keyword or identifier, found "
        << currToken;
      }
    catch(vtkFoamError& e)
      {
      if(isSubDictionary)
        {
        throw;
        }
      else
        {
        io.setError(e);
        return false;
        }
      }
  }
};

void vtkOpenFOAMReaderPrivate::vtkFoamIOobject::readHeader()
{
  vtkFoamToken firstToken;

  readExpecting("FoamFile");
  readExpecting('{');

  vtkFoamDict headerDict;
  // throw exception in case of error
  headerDict.read(*this, true, vtkFoamToken());

  const vtkFoamEntry& formatEntry = headerDict.lookup("format");
  if(!formatEntry.found())
    {
    throw vtkFoamError()
      << "format entry (binary/ascii) not found in FoamFile header";
    }
  // case does matter (e. g. "BINARY" is treated as ascii)
  // cf. src/OpenFOAM/db/IOstreams/IOstreams/IOstream.C
  format_ = (formatEntry.toString() == "binary" ? BINARY : ASCII);

  const vtkFoamEntry& classEntry = headerDict.lookup("class");
  if(!classEntry.found())
    {
    throw vtkFoamError() << "class name not found in FoamFile header";
    }
  headerClassName_ = classEntry.toString();

  const vtkFoamEntry& objectEntry = headerDict.lookup("object");
  if(!objectEntry.found())
    {
    throw vtkFoamError() << "object name not found in FoamFile header";
    }
  objectName_ = objectEntry.toString();
}

vtkOpenFOAMReaderPrivate::vtkFoamEntryValue::vtkFoamEntryValue(
  vtkFoamEntryValue& value, vtkFoamEntry *upperEntryPtr)
  : vtkFoamToken(value), isUniform_(value.isUniform()), managed_(true),
    upperEntryPtr_(upperEntryPtr)
{
  switch(type_)
    {
    case LABELLIST: case SCALARLIST: case VECTORLIST: case STRINGLIST:
      vtkObjectPtr_ = value.toVTKObject();
      vtkObjectPtr_->Register(0); // increase reference count +1
      break;
    case LABELLISTLIST:
      labelListListPtr_ = new intVectorVector(*value.labelListListPtr_);
      break;
    case ENTRYVALUELIST:
      {
      const int nValues = value.entryValuePtrs_->size();
      entryValuePtrs_ = new vtkstd::vector<vtkFoamEntryValue*>(nValues);
      for(int valueI = 0; valueI < nValues; valueI++)
	{
        entryValuePtrs_->operator[](valueI) = new vtkFoamEntryValue(
          *value.entryValuePtrs_->operator[](valueI), upperEntryPtr_);
	}
      }
      break;
    case DICTIONARY:
      // upperEntryPtr_ is null when called from vtkFoamDict constructor
      if(upperEntryPtr_ != NULL)
        {
        dictPtr_
          = new vtkFoamDict(*value.dictPtr_, upperEntryPtr_->upperDictPtr());
        }
      else
        {
        dictPtr_ = NULL;
        }
      break;
    case EMPTYLIST:
      break;
    }
}

void vtkOpenFOAMReaderPrivate::vtkFoamEntryValue::clear()
{
  if(managed_)
    {
    switch(type_)
      {
      case LABELLIST: case SCALARLIST: case VECTORLIST: case STRINGLIST:
        vtkObjectPtr_->Delete();
        break;
      case LABELLISTLIST:
        delete labelListListPtr_;
        break;
      case ENTRYVALUELIST:
        for(size_t valueI = 0; valueI < entryValuePtrs_->size() ; valueI++)
          {
          delete entryValuePtrs_->operator[](valueI);
          }
        delete entryValuePtrs_;
        break;
      case DICTIONARY:
        delete dictPtr_;
        break;
      }
    }
}

// general-purpose list reader - guess the type of the list and read
// it. only supports ascii format and assumes the preceding '(' has
// already been thrown away.  the reader supports nested list with
// variable lengths (e. g. `((token token) (token token token)).'
// also supports compound of tokens and lists (e. g. `((token token)
// token)') only if a list comes as the first value.
void vtkOpenFOAMReaderPrivate::vtkFoamEntryValue::readList(vtkFoamIOobject& io)
{
  vtkFoamToken currToken;
  io.read(currToken);

  // initial guess of the list type
  if(currToken.type() == vtkFoamToken::LABEL)
    {
    // if the first token is of type LABEL it might be either an element of
    // a labelList or the size of a sublist so proceed to the next token
    vtkFoamToken nextToken;
    if(!io.read(nextToken))
      {
      throw vtkFoamError() << "Unexpected EOF";
      }
    if(nextToken.type() == vtkFoamToken::LABEL)
      {
      labelListPtr_ = vtkIntArray::New();
      labelListPtr_->InsertNextValue(currToken.to<int>());
      labelListPtr_->InsertNextValue(nextToken.to<int>());
      type_ = LABELLIST;
      }
    else if(nextToken.type() == vtkFoamToken::SCALAR)
      {
      scalarListPtr_ = vtkFloatArray::New();
      scalarListPtr_->InsertNextValue(currToken.to<float>());
      scalarListPtr_->InsertNextValue(nextToken.to<float>());
      type_ = SCALARLIST;
      }
    else if(nextToken == '(') // list of list: read recursively
      {
      entryValuePtrs_ = new vtkstd::vector<vtkFoamEntryValue*>;
      entryValuePtrs_->push_back(new vtkFoamEntryValue(upperEntryPtr_));
      entryValuePtrs_->back()->readList(io);
      type_ = ENTRYVALUELIST;
      }
    else if(nextToken == ')') // list with only one label element
      {
      labelListPtr_ = vtkIntArray::New();
      labelListPtr_->SetNumberOfValues(1);
      labelListPtr_->SetValue(0, currToken.to<int>());
      type_ = LABELLIST;
      return;
      }
    else
      {
      throw vtkFoamError() << "Expected number, '(' or ')', found "
        << nextToken;
      }
    }
  else if(currToken.type() == vtkFoamToken::SCALAR)
    {
    scalarListPtr_ = vtkFloatArray::New();
    scalarListPtr_->InsertNextValue(currToken.to<float>());
    type_ = SCALARLIST;
    }
  // if the first word is a string we have to read another token to determine
  // if the first word is a keyword for the following dictionary
  else if(currToken.type() == vtkFoamToken::STRING)
    {
    vtkFoamToken nextToken;
    if(!io.read(nextToken))
      {
      throw vtkFoamError() << "Unexpected EOF";
      }
    if(nextToken.type() == vtkFoamToken::STRING) // list of strings
      {
      stringListPtr_ = vtkStringArray::New();
      stringListPtr_->InsertNextValue(currToken.toString());
      stringListPtr_->InsertNextValue(nextToken.toString());
      type_ = STRINGLIST;
      }
    // dictionary with the already read stringToken as the first keyword
    else if(nextToken == '{')
      {
      if(currToken.toString() == "")
        {
        throw "Empty string is invalid as a keyword for dictionary entry";
        }
      readDictionary(io, currToken);
      // the dictionary read as list has the entry terminator ';' so
      // we have to skip it
      return;
      }
    else if(nextToken == ')') // list with only one string element
      {
      stringListPtr_ = vtkStringArray::New();
      stringListPtr_->SetNumberOfValues(1);
      stringListPtr_->SetValue(0, currToken.toString());
      type_ = STRINGLIST;
      return;
      }
    else
      {
      throw vtkFoamError() << "Expected string, '{' or ')', found "
        << nextToken;
      }
    }
  else if(currToken == '(') // list of lists: read recursively
    {
    entryValuePtrs_ = new vtkstd::vector<vtkFoamEntryValue*>;
    entryValuePtrs_->push_back(new vtkFoamEntryValue(upperEntryPtr_));
    entryValuePtrs_->back()->readList(io);
    // read all the following values as arbitrary entryValues
    // the alphaContactAngle b.c. in multiphaseInterFoam/damBreak4phase
    // reaquires this treatment (reading by readList() is not enough)
    do
      {
      entryValuePtrs_->push_back(new vtkFoamEntryValue(upperEntryPtr_));
      entryValuePtrs_->back()->read(io);
      }
    while(*entryValuePtrs_->back() != ')' && *entryValuePtrs_->back() != '}'
      && *entryValuePtrs_->back() != ';');

    if(*entryValuePtrs_->back() != ')')
      {
      throw vtkFoamError() << "Expected ')' before "
        << *entryValuePtrs_->back();
      }

    // delete ')'
    delete entryValuePtrs_->back();
    entryValuePtrs_->pop_back();
    type_ = ENTRYVALUELIST;
    return;
    }
  else if(currToken == ')') // empty list
    {
    type_ = EMPTYLIST;
    return;
    }
  // FIXME: may (or may not) need identifier handling

  while(io.read(currToken) && currToken != ')')
    {
    if(type_ == LABELLIST)
      {
      if(currToken.type() == vtkFoamToken::SCALAR) // switch to scalarList
        {
        // labelListPtr_ and scalarListPtr_ are packed into a single union so
        // we need a temprary pointer
        vtkFloatArray* slPtr = vtkFloatArray::New();
        const int size = labelListPtr_->GetNumberOfTuples();
        slPtr->SetNumberOfValues(size + 1);
        for(int i = 0; i < size; i++)
          {
          slPtr->SetValue(i, static_cast<float>(labelListPtr_->GetValue(i)));
          }
        labelListPtr_->Delete();
        slPtr->SetValue(size, currToken.to<float>());
        scalarListPtr_ = slPtr; // copy after labelListPtr_ is deleted
        type_ = SCALARLIST;
        }
      else if(currToken.type() == vtkFoamToken::LABEL)
        {
	labelListPtr_->InsertNextValue(currToken.to<int>());
        }
      else
        {
        throw vtkFoamError() << "Expected a number, found " << currToken;
        }
      }
    else if(type_ == SCALARLIST)
      {
      if(currToken.is<float>())
        {
	scalarListPtr_->InsertNextValue(currToken.to<float>());
        }
      else
        {
        throw vtkFoamError() << "Expected a number, found " << currToken;
        }
      }
    else if(type_ == STRINGLIST)
      {
      if(currToken.type() == vtkFoamToken::STRING)
        {
        stringListPtr_->InsertNextValue(currToken.toString());
        }
      else
        {
        throw vtkFoamError() << "Expected a string, found " << currToken;
        }
      }
    else if(type_ == ENTRYVALUELIST)
      {
      if(currToken.type() == vtkFoamToken::LABEL)
        {
        // skip the number of elements to make things simple
        if(!io.read(currToken))
          {
          throw vtkFoamError() << "Unexpected EOF";
          }
        }
      if(currToken != '(')
        {
        throw vtkFoamError() << "Expected '(', found " << currToken;
        }
      entryValuePtrs_->push_back(new vtkFoamEntryValue(upperEntryPtr_));
      entryValuePtrs_->back()->readList(io);
      }
    else
      {
      throw vtkFoamError() << "Unexpected token " << currToken;
      }
    }

  if(type_ == LABELLIST)
    {
    labelListPtr_->Squeeze();
    }
  else if(type_ == SCALARLIST)
    {
    scalarListPtr_->Squeeze();
    }
  else if(type_ == STRINGLIST)
    {
    stringListPtr_->Squeeze();
    }
}

// a list of dictionaries is actually read as a dictionary
void vtkOpenFOAMReaderPrivate::vtkFoamEntryValue::readDictionary(
  vtkFoamIOobject& io, const vtkFoamToken& firstKeyword)
{
  dictPtr_ = new vtkFoamDict(upperEntryPtr_->upperDictPtr());
  type_ = DICTIONARY;
  dictPtr_->read(io, true, firstKeyword);
}

// guess the type of the given entry value and read it
void vtkOpenFOAMReaderPrivate::vtkFoamEntryValue::read(vtkFoamIOobject& io)
{
  vtkFoamToken currToken;
  if(!io.read(currToken))
    {
    throw vtkFoamError() << "Unexpected EOF";
    }

  if(currToken == '{')
    {
    readDictionary(io, vtkFoamToken());
    return;
    }
  // for reading sublist from vtkFoamEntryValue::readList() or there
  // are cases where lists without the (non)uniform keyword appear
  // (e. g. coodles/pitsDaily/0/U, uniformFixedValue b.c.)
  else if(currToken == '(')
    {
    readList(io);
    return;
    }
  else if(currToken == '[')
    {
    readDimensionSet(io);
    return;
    }
  else if(currToken == "uniform")
    {
    if(!io.read(currToken))
      {
      throw vtkFoamError()
        << "Expected a uniform value or a list, found unexpected EOF";
      }
    if(currToken == '(')
      {
      readList(io);
      }
    else if(currToken.type() == vtkFoamToken::LABEL
      || currToken.type() == vtkFoamToken::SCALAR
      || currToken.type() == vtkFoamToken::STRING)
      {
      vtkFoamToken::operator=(currToken);
      }
    else // unexpected punctuation token
      {
      throw vtkFoamError() << "Expected number, string or (, found "
        << currToken;
      }
    isUniform_ = true;
    }
  else if(currToken == "nonuniform")
    {
    if(!io.read(currToken))
      {
      throw vtkFoamError() << "Expected list type specifier, found EOF";
      }
    isUniform_ = false;
    if(currToken == "List<scalar>")
      {
      readNonuniformList<SCALARLIST, listTraits<vtkFloatArray, float> >(io);
      }
    else if(currToken == "List<sphericalTensor>")
      {
      readNonuniformList<VECTORLIST,
        vectorListTraits<vtkFloatArray, float, 1, false> >(io);
      }
    else if(currToken == "List<vector>")
      {
      readNonuniformList<VECTORLIST,
        vectorListTraits<vtkFloatArray, float, 3, false> >(io);
      }
    else if(currToken == "List<symmTensor>")
      {
      readNonuniformList<VECTORLIST,
        vectorListTraits<vtkFloatArray, float, 6, false> >(io);
      }
    else if(currToken == "List<tensor>")
      {
      readNonuniformList<VECTORLIST,
	vectorListTraits<vtkFloatArray, float, 9, false> >(io);
      }
    // List<bool> is read as List<label>
    else if(currToken =="List<label>" || currToken == "List<bool>")
      {
      readNonuniformList<LABELLIST, listTraits<vtkIntArray, int> >(io);
      }
    // an empty list doesn't have a list type specifier
    else if(currToken.type() == vtkFoamToken::LABEL
      && currToken.to<int>() == 0)
      {
      type_ = EMPTYLIST;
      io.readExpecting('(');
      io.readExpecting(')');
      }
    else
      {
      throw vtkFoamError() << "Unsupported nonuniform list type " << currToken;
      }
    }
  // zones have list without a uniform/nonuniform keyword
  // List<bool> is read as List<label>
  // (e. g. flipMap entry in faceZones)
  else if(currToken == "List<label>" || currToken == "List<bool>")
    {
    isUniform_ = false;
    readNonuniformList<LABELLIST, listTraits<vtkIntArray, int> >(io);
    }
  else if(currToken.type() == vtkFoamToken::PUNCTUATION
    || currToken.type() == vtkFoamToken::LABEL
    || currToken.type() == vtkFoamToken::SCALAR
    || currToken.type() == vtkFoamToken::STRING
    || currToken.type() == vtkFoamToken::IDENTIFIER)
    {
    vtkFoamToken::operator=(currToken);
    }
}

// read values of an entry
void vtkOpenFOAMReaderPrivate::vtkFoamEntry::read(vtkFoamIOobject& io)
{
  for(;;)
    {
    Superclass::push_back(new vtkFoamEntryValue(this));
    Superclass::back()->read(io);

    if(Superclass::size() >= 2)
      {
      vtkFoamEntryValue& secondLastValue
        = *Superclass::operator[](Superclass::size() - 2);
      if(secondLastValue.type() == vtkFoamToken::LABEL)
        {
        vtkFoamEntryValue& lastValue = *Superclass::back();

        // a zero-sized nonuniform list without prefixing "nonuniform"
        // keyword nor list type specifier (i. e. `0()';
        // e. g. simpleEngine/0/polyMesh/pointZones) requires special
        // care (one with nonuniform prefix is treated within
        // vtkFoamEntryValue::read()). still this causes errornous
        // behavior for `0 nonuniform 0()' but this should be extremely
        // rare
        if(lastValue.type() == vtkFoamToken::EMPTYLIST
          && secondLastValue == 0)
          {
          delete Superclass::back();
          Superclass::pop_back(); // delete the last value
          Superclass::back()->setEmptyList(); // mark new last value as empty
          }
        // for an exceptional expression of `LABEL{LABELorSCALAR}' without
        // type prefix (e. g. `2{-0}' in mixedRhoE B.C. in
        // rhopSonicFoam/shockTube)
        else if(lastValue.type() == vtkFoamToken::DICTIONARY)
          {
          if(lastValue.dictionary().type() == vtkFoamToken::UNIFORMLABELLIST)
            {
            const int size = secondLastValue.to<int>();
            const int value = lastValue.dictionary().to<int>();
            // delete last two values
            delete Superclass::back();
            Superclass::pop_back();
            delete Superclass::back();
            Superclass::pop_back();
            // make new labelList
            Superclass::push_back(new vtkFoamEntryValue(this));
            Superclass::back()->makeLabelList(value, size);
            }
          else if(lastValue.dictionary().type()
            == vtkFoamToken::UNIFORMSCALARLIST)
            {
	    const int size = secondLastValue.to<int>();
            const float value = lastValue.dictionary().to<float>();
            // delete last two values
            delete Superclass::back();
            Superclass::pop_back();
            delete Superclass::back();
            Superclass::pop_back();
            // make new labelList
            Superclass::push_back(new vtkFoamEntryValue(this));
            Superclass::back()->makeScalarList(value, size);
            }
          }
        }
      }

    if(Superclass::back()->type() == vtkFoamToken::IDENTIFIER)
      {
      // substitute identifier
      const vtkStdString identifier(Superclass::back()->toIdentifier());
      delete Superclass::back();
      Superclass::pop_back();

      for(vtkFoamDict *uDictPtr = upperDictPtr_;;)
        {
        const vtkFoamEntry& identifiedEntry = uDictPtr->lookup(identifier);

        if(identifiedEntry.found())
          {
          for(size_t valueI = 0; valueI < identifiedEntry.size(); valueI++)
            {
            Superclass::push_back(
              new vtkFoamEntryValue(identifiedEntry.value(valueI), this));
            }
          break;
          }
        else
          {
          uDictPtr = uDictPtr->upperDictPtr();
          if(uDictPtr == NULL)
            {
            throw vtkFoamError() << "substituting entry " << identifier
		 << " not found";
            }
          }
        }
      }
    else if(*Superclass::back() == ';')
      {
      delete Superclass::back();
      Superclass::pop_back();
      break;
      }
    else if(Superclass::back()->type() == vtkFoamToken::DICTIONARY)
      {
      // subdictionary is not suffixed by an entry terminator ';'
      break;
      }
    else if(*Superclass::back() == '}' || *Superclass::back() == ')')
      {
      throw vtkFoamError() << "Unmatched " << *Superclass::back();
      }
    }
}

//-----------------------------------------------------------------------------
// vtkOpenFOAMReaderPrivate constructor and destructor
vtkOpenFOAMReaderPrivate::vtkOpenFOAMReaderPrivate()
{
  // DATA TIMES
  this->TimeStep = 0;
  this->TimeStepOld = -1;
  this->TimeValues = vtkDoubleArray::New();
  this->TimeNames = vtkStringArray::New();

  // selection
  this->InternalMeshSelectionStatus = 0;
  this->InternalMeshSelectionStatusOld = 0;

  // DATA COUNTS
  this->NumCells = 0;
  this->NumPoints = 0;

  this->VolFieldFiles = vtkStringArray::New();
  this->PointFieldFiles = vtkStringArray::New();
  this->LagrangianFieldFiles = vtkStringArray::New();
  this->PolyMeshPointsDir = vtkStringArray::New();
  this->PolyMeshFacesDir = vtkStringArray::New();

  // for creating cell-to-point translated data
  this->BoundaryPointMap = NULL;
  this->AllBoundaries = NULL;
  this->AllBoundariesPointMap = NULL;
  this->InternalPoints = NULL;

  // for caching mesh
  this->InternalMesh = NULL;
  this->BoundaryMesh = NULL;
  this->BoundaryPointMap = NULL;
  this->FaceOwner = NULL;
  this->PointZoneMesh = NULL;
  this->FaceZoneMesh = NULL;
  this->CellZoneMesh = NULL;

  // for decomposing polyhedra
  this->AdditionalCellIds = NULL;
  this->AdditionalCellPoints = NULL;
}

vtkOpenFOAMReaderPrivate::~vtkOpenFOAMReaderPrivate()
{
  this->TimeValues->Delete();
  this->TimeNames->Delete();

  this->PolyMeshPointsDir->Delete();
  this->PolyMeshFacesDir->Delete();
  this->VolFieldFiles->Delete();
  this->PointFieldFiles->Delete();
  this->LagrangianFieldFiles->Delete();

  this->ClearMeshes();
}

void vtkOpenFOAMReaderPrivate::ClearInternalMeshes()
{
  if(this->FaceOwner != NULL)
    {
    this->FaceOwner->Delete();
    this->FaceOwner = NULL;
    }
  if(this->InternalMesh != NULL)
    {
    this->InternalMesh->Delete();
    this->InternalMesh = NULL;
    }
  if(this->AdditionalCellIds != NULL)
    {
    this->AdditionalCellIds->Delete();
    this->AdditionalCellIds = NULL;
    }
  delete this->AdditionalCellPoints;
  this->AdditionalCellPoints = NULL;

  if(this->PointZoneMesh != NULL)
    {
    this->PointZoneMesh->Delete();
    this->PointZoneMesh = NULL;
    }
  if(this->FaceZoneMesh != NULL)
    {
    this->FaceZoneMesh->Delete();
    this->FaceZoneMesh = NULL;
    }
  if(this->CellZoneMesh != NULL)
    {
    this->CellZoneMesh->Delete();
    this->CellZoneMesh = NULL;
    }
}

void vtkOpenFOAMReaderPrivate::ClearBoundaryMeshes()
{
  if(this->BoundaryMesh != NULL)
    {
    this->BoundaryMesh->Delete();
    this->BoundaryMesh = NULL;
    }

  delete this->BoundaryPointMap;
  this->BoundaryPointMap = NULL;

  if(this->InternalPoints != NULL)
    {
    this->InternalPoints->Delete();
    this->InternalPoints = NULL;
    }
  if(this->AllBoundaries != NULL)
    {
    this->AllBoundaries->Delete();
    this->AllBoundaries = NULL;
    }
  if(this->AllBoundariesPointMap != NULL)
    {
    this->AllBoundariesPointMap->Delete();
    this->AllBoundariesPointMap = NULL;
    }
}

void vtkOpenFOAMReaderPrivate::ClearMeshes()
{
  this->ClearInternalMeshes();
  this->ClearBoundaryMeshes();
}

void vtkOpenFOAMReaderPrivate::SetTimeValue(const double requestedTime)
{
  const int nTimes = this->TimeValues->GetNumberOfTuples();
  if(nTimes > 0)
    {
    int minTimeI = 0;
    double minTimeDiff = fabs(this->TimeValues->GetValue(0) - requestedTime);
    for(int timeI = 1; timeI < nTimes; timeI++)
      {
      const double timeDiff(
        fabs(this->TimeValues->GetValue(timeI) - requestedTime));
      if(timeDiff < minTimeDiff)
        {
        minTimeI = timeI;
        minTimeDiff = timeDiff;
        }
      }
    this->SetTimeStep(minTimeI); // set Modified() if TimeStep changed
    }
}

//-----------------------------------------------------------------------------
void vtkOpenFOAMReaderPrivate::SetupInformation(
  const vtkStdString &casePath, const vtkStdString &regionName,
  const vtkStdString &procName, vtkOpenFOAMReaderPrivate *master)
{
  // copy parent, path and timestep information from master
  this->CasePath = casePath;
  this->RegionName = regionName;
  this->ProcessorName = procName;
  this->Parent = master->Parent;
  this->TimeValues->Delete();
  this->TimeValues = master->TimeValues;
  this->TimeValues->Register(0);
  this->TimeNames->Delete();
  this->TimeNames = master->TimeNames;
  this->TimeNames->Register(0);

  this->PopulatePolyMeshDirArrays();
}

//-----------------------------------------------------------------------------
void vtkOpenFOAMReaderPrivate::GetFieldNames(const vtkStdString &tempPath,
  const bool isLagrangian, vtkStringArray *cellObjectNames,
  vtkStringArray *pointObjectNames)
{
  // open the directory and get num of files
  vtkDirectory *directory = vtkDirectory::New();
  if(!directory->Open(tempPath.c_str()))
    {
    // no data
    directory->Delete();
    return;
    }

  // loop over all files and locate valid fields
  int nFieldFiles = directory->GetNumberOfFiles();
  for(int j = 0; j < nFieldFiles; j++)
    {
    const vtkStdString fieldFile(directory->GetFile(j));
    const int len = fieldFile.length();

    // excluded extensions cf. src/OpenFOAM/OSspecific/Unix/Unix.C
    if(!directory->FileIsDirectory(fieldFile.c_str())
      && fieldFile.substr(len - 1) != "~"
      && (len < 4 || (fieldFile.substr(len - 4) != ".bak"
      && fieldFile.substr(len - 4) != ".BAK"
      && fieldFile.substr(len - 4) != ".old"))
      && (len < 5 || fieldFile.substr(len - 5) != ".save"))
      {
      vtkFoamIOobject io(this->CasePath);
      if(io.open(tempPath + "/" + fieldFile)) // file exists and readable
        {
        const vtkStdString& cn = io.className();
        if(isLagrangian)
          {
          if(cn == "scalarField" || cn == "vectorField"
            || cn == "sphericalTensorField" || cn == "symmTensorField"
            || cn == "tensorField")
            {
            // real file name
            this->LagrangianFieldFiles->InsertNextValue(fieldFile);
            // object name
            pointObjectNames->InsertNextValue(io.objectName());
            }
          }
        else
          {
          if(cn == "volScalarField" || cn == "pointScalarField"
            || cn == "volVectorField" || cn == "pointVectorField"
            || cn == "volSphericalTensorField"
            || cn == "pointSphericalTensorField"
            || cn == "volSymmTensorField" || cn == "pointSymmTensorField"
            || cn == "volTensorField" || cn == "pointTensorField")
            {
            if(cn.substr(0, 3) == "vol")
              {
              // real file name
              this->VolFieldFiles->InsertNextValue(fieldFile);
              // object name
              cellObjectNames->InsertNextValue(io.objectName());
              }
            else
              {
              this->PointFieldFiles->InsertNextValue(fieldFile);
              pointObjectNames->InsertNextValue(io.objectName());
              }
            }
          }
        io.close();
        }
      }
    }
  // inserted objects are squeezed later in SortFieldFiles()
  directory->Delete();
}

//-----------------------------------------------------------------------------
// locate laglangian clouds
void vtkOpenFOAMReaderPrivate::LocateLagrangianClouds(
  vtkStringArray *lagrangianObjectNames, const vtkStdString &timePath)
{
  vtkDirectory *directory = vtkDirectory::New();
  if(directory->Open((timePath + this->RegionPath() + "/lagrangian").c_str()))
    {
    // search for sub-clouds (OF 1.5 format)
    const int nFiles = directory->GetNumberOfFiles();
    bool isSubCloud = false;
    for(int fileI = 0; fileI < nFiles; fileI++)
      {
      const vtkStdString fileNameI(directory->GetFile(fileI));
      if(fileNameI != "." && fileNameI != ".."
	&& directory->FileIsDirectory(fileNameI.c_str()))
        {
        vtkFoamIOobject io(this->CasePath);
        const vtkStdString subCloudName(this->RegionPrefix() + "lagrangian/"
          + fileNameI);
        const vtkStdString subCloudFullPath(timePath + "/" + subCloudName);
        // lagrangian positions. there are many concrete class names
        // e. g. Cloud<parcel>, basicKinematicCloud etc.
        if(io.open(subCloudFullPath + "/positions")
          && io.className().find("Cloud") != vtkStdString::npos
          && io.objectName() == "positions")
          {
          isSubCloud = true;
	  // a lagrangianPath has to be in a bit different format from
	  // subCloudName to make the "lagrangian" reserved path
	  // component and a mesh region with the same name
	  // distinguishable later
	  const vtkStdString subCloudPath(
	    this->RegionName + "/lagrangian/" + fileNameI);
          if(this->Parent->LagrangianPaths->LookupValue(subCloudPath) == -1)
            {
            this->Parent->LagrangianPaths->InsertNextValue(subCloudPath);
            }
          this->GetFieldNames(subCloudFullPath, true, NULL,
            lagrangianObjectNames);
	  this->Parent->PatchDataArraySelection->AddArray(subCloudName.c_str());
          }
        }
      }
    // if there's no sub-cloud then OF < 1.5 format
    if(!isSubCloud)
      {
      vtkFoamIOobject io(this->CasePath);
      const vtkStdString cloudName(this->RegionPrefix() + "lagrangian");
      const vtkStdString cloudFullPath(timePath + "/" + cloudName);
      if(io.open(cloudFullPath + "/positions")
        && io.className() == "Cloud" && io.objectName() == "positions")
        {
        const vtkStdString cloudPath(this->RegionName + "/lagrangian");
        if(this->Parent->LagrangianPaths->LookupValue(cloudPath) == -1)
          {
          this->Parent->LagrangianPaths->InsertNextValue(cloudPath);
          }
        this->GetFieldNames(cloudFullPath, true, NULL, lagrangianObjectNames);
	this->Parent->PatchDataArraySelection->AddArray(cloudName.c_str());
        }
      }
    this->Parent->LagrangianPaths->Squeeze();
    }
  directory->Delete();
}

//-----------------------------------------------------------------------------
void vtkOpenFOAMReaderPrivate::SortFieldFiles(vtkStringArray *selections,
  vtkStringArray *files, vtkStringArray *objects)
{
  objects->Squeeze();
  files->Squeeze();
  vtkSortDataArray::Sort(objects, files);
  for(int nameI = 0; nameI < objects->GetNumberOfValues(); nameI++)
    {
    selections->InsertNextValue(objects->GetValue(nameI));
    }
  objects->Delete();
}

//-----------------------------------------------------------------------------
// create field data lists and cell/point array selection lists
int vtkOpenFOAMReaderPrivate::MakeMetaDataAtTimeStep(
  vtkStringArray *cellSelectionNames, vtkStringArray *pointSelectionNames,
  vtkStringArray *lagrangianSelectionNames, const bool listNextTimeStep)
{
  // Read the patches from the boundary file into selection array
  if(this->PolyMeshFacesDir->GetValue(this->TimeStep)
    != this->BoundaryDict.timeDir
    || this->Parent->PatchDataArraySelection->GetMTime()
    != this->Parent->PatchSelectionMTimeOld)
    {
    this->BoundaryDict.clear();
    this->BoundaryDict.timeDir
      = this->PolyMeshFacesDir->GetValue(this->TimeStep);

    const bool isSubRegion = this->RegionName != "";
    vtkFoamDict *boundaryDict = this->GatherBlocks("boundary", isSubRegion);
    if(boundaryDict == NULL)
      {
      if(isSubRegion)
        {
        return 0;
        }
      }
    else
      {
      // Add the internal mesh by default always
      const vtkStdString internalMeshName(
        this->RegionPrefix() + "internalMesh");
      this->Parent->PatchDataArraySelection->AddArray(internalMeshName.c_str());
      this->InternalMeshSelectionStatus
        = this->Parent->GetPatchArrayStatus(internalMeshName.c_str());

      // iterate through each entry in the boundary file
      int allBoundariesNextStartFace = 0;
      this->BoundaryDict.resize(boundaryDict->size());
      for(size_t i = 0; i < boundaryDict->size(); i++)
        {
        vtkFoamEntry &boundaryEntryI = boundaryDict->entry(i);
        vtkFoamEntry &nFacesEntry
          = boundaryEntryI.dictionary().lookup("nFaces");
        if(!nFacesEntry.found())
          {
          vtkErrorMacro(<< "nFaces entry not found in boundary entry "
            << boundaryEntryI.keyword().c_str());
          delete boundaryDict;
          return 0;
          }
        const int nFaces = nFacesEntry.toInt();

        // extract name of the current patch for insertion
        const vtkStdString &boundaryNameI = boundaryEntryI.keyword();

        // create BoundaryDict entry
        vtkFoamBoundaryEntry &BoundaryDictRef = this->BoundaryDict[i];
        BoundaryDictRef.nFaces = nFaces;
        BoundaryDictRef.boundaryName = boundaryNameI;
        vtkFoamEntry &startFaceEntry
          = boundaryEntryI.dictionary().lookup("startFace");
        if(!startFaceEntry.found())
          {
          vtkErrorMacro(<< "startFace entry not found in boundary entry "
            << boundaryEntryI.keyword().c_str());
          delete boundaryDict;
          return 0;
          }
        BoundaryDictRef.startFace = startFaceEntry.toInt();
        vtkFoamEntry &typeEntry = boundaryEntryI.dictionary().lookup("type");
        if(!typeEntry.found())
          {
          vtkErrorMacro(<< "type entry not found in boundary entry "
            << boundaryEntryI.keyword().c_str());
          delete boundaryDict;
          return 0;
          }
        BoundaryDictRef.allBoundariesStartFace = allBoundariesNextStartFace;
        const vtkStdString &typeNameI = typeEntry.toString();
        // if the basic type of the patch is one of the followings the
        // point-filtered values at patches are overrided by patch values
        if(typeNameI == "patch" || typeNameI == "wall")
          {
          BoundaryDictRef.boundaryType = vtkFoamBoundaryEntry::PHYSICAL;
          allBoundariesNextStartFace += nFaces;
          }
        else if(typeNameI == "processor")
          {
          BoundaryDictRef.boundaryType = vtkFoamBoundaryEntry::PROCESSOR;
          allBoundariesNextStartFace += nFaces;
          }
        else
          {
          BoundaryDictRef.boundaryType = vtkFoamBoundaryEntry::GEOMETRICAL;
          }
        BoundaryDictRef.isActive = false;

        // always hide processor patches for decomposed cases to keep
        // vtkAppendCompositeDataLeaves happy
        if(this->ProcessorName != ""
          && BoundaryDictRef.boundaryType == vtkFoamBoundaryEntry::PROCESSOR)
          {
          continue;
          }
        const vtkStdString selectionName(this->RegionPrefix() + boundaryNameI);
        if(this->Parent->PatchDataArraySelection->
          ArrayExists(selectionName.c_str()))
          {
          // Mark boundary if selected for display
          if(this->Parent->GetPatchArrayStatus(selectionName.c_str()))
            {
            BoundaryDictRef.isActive = true;
            }
          }
        else
          {
          // add patch to list with selection status turned off:
          // the patch is added to list even if its size is zero
          this->Parent->PatchDataArraySelection->DisableArray(
            selectionName.c_str());
          }
        }

      delete boundaryDict;
      }
    }

  // Add scalars and vectors to metadata
  vtkStdString timePath(this->CurrentTimePath());
  // do not do "RemoveAllArrays()" to accumulate array selections
  // this->CellDataArraySelection->RemoveAllArrays();
  this->VolFieldFiles->Initialize();
  this->PointFieldFiles->Initialize();
  vtkStringArray *cellObjectNames = vtkStringArray::New();
  vtkStringArray *pointObjectNames = vtkStringArray::New();
  this->GetFieldNames(timePath + this->RegionPath(), false, cellObjectNames,
    pointObjectNames);

  this->LagrangianFieldFiles->Initialize();
  if(listNextTimeStep)
    {
    this->Parent->LagrangianPaths->Initialize();
    }
  vtkStringArray *lagrangianObjectNames = vtkStringArray::New();
  this->LocateLagrangianClouds(lagrangianObjectNames, timePath);

  // if the requested timestep is 0 then we also look at the next
  // timestep to add extra objects that don't exist at timestep 0 into
  // selection lists. Note the ObjectNames array will be recreated in
  // RequestData() so we don't have to worry about duplicated fields.
  if(listNextTimeStep && this->TimeValues->GetNumberOfTuples() >= 2
    && this->TimeStep == 0)
    {
    const vtkStdString timePath(this->TimePath(1));
    this->GetFieldNames(timePath + this->RegionPath(), false, cellObjectNames,
      pointObjectNames);
    // if lagrangian clouds were not found at timestep 0
    if(this->Parent->LagrangianPaths->GetNumberOfTuples() == 0)
      {
      this->LocateLagrangianClouds(lagrangianObjectNames, timePath);
      }
    }

  // sort array names
  this->SortFieldFiles(cellSelectionNames, this->VolFieldFiles,
    cellObjectNames);
  this->SortFieldFiles(pointSelectionNames, this->PointFieldFiles,
    pointObjectNames);
  this->SortFieldFiles(lagrangianSelectionNames, this->LagrangianFieldFiles,
    lagrangianObjectNames);

  return 1;
}

//-----------------------------------------------------------------------------
// list time directories according to controlDict
bool vtkOpenFOAMReaderPrivate::ListTimeDirectoriesByControlDict(
  vtkFoamDict* dictPtr)
{
  vtkFoamDict& dict = *dictPtr;

  vtkFoamEntry& startTimeEntry = dict.lookup("startTime");
  if(!startTimeEntry.found())
    {
    vtkErrorMacro(<< "startTime entry not found in controlDict");
    return false;
    }
  // using double to precisely handle time values
  const double startTime = startTimeEntry.toDouble();

  vtkFoamEntry& endTimeEntry = dict.lookup("endTime");
  if(!endTimeEntry.found())
    {
    vtkErrorMacro(<< "endTime entry not found in controlDict");
    return false;
    }
  const double endTime = endTimeEntry.toDouble();

  vtkFoamEntry& deltaTEntry = dict.lookup("deltaT");
  if(!deltaTEntry.found())
    {
    vtkErrorMacro(<< "deltaT entry not found in controlDict");
    return false;
    }
  const double deltaT = deltaTEntry.toDouble();

  vtkFoamEntry& writeIntervalEntry = dict.lookup("writeInterval");
  if(!writeIntervalEntry.found())
    {
    vtkErrorMacro(<< "writeInterval entry not found in controlDict");
    return false;
    }
  const double writeInterval = writeIntervalEntry.toDouble();

  vtkFoamEntry& timeFormatEntry = dict.lookup("timeFormat");
  if(!timeFormatEntry.found())
    {
    vtkErrorMacro(<< "timeFormat entry not found in controlDict");
    return false;
    }
  const vtkStdString timeFormat(timeFormatEntry.toString());

  vtkFoamEntry& timePrecisionEntry = dict.lookup("timePrecision");
  const int timePrecision  // default is 6
    = (timePrecisionEntry.found() ? timePrecisionEntry.toInt() : 6);

  // calculate the time step increment based on type of run
  const vtkStdString writeControl(dict.lookup("writeControl").toString());
  double timeStepIncrement;
  if(writeControl == "timeStep")
    {
    timeStepIncrement = writeInterval * deltaT;
    }
  else if(writeControl == "runTime" || writeControl == "adjustableRunTime")
    {
    timeStepIncrement = writeInterval;
    }
  else
    {
    vtkErrorMacro(<<"Time step can't be determined because writeControl is"
      " set to " << writeControl.c_str());
    return false;
    }

  // calculate how many timesteps there should be
  const double tempResult = (endTime - startTime) / timeStepIncrement;
  // +0.5 to round up
  const int tempNumTimeSteps = static_cast<int>(tempResult + 0.5) + 1;

  // make sure time step dir exists
  vtkstd::vector<double> tempSteps;
  vtkDirectory *test = vtkDirectory::New();
  this->TimeValues->Initialize();
  this->TimeNames->Initialize();

  // determine time name based on Foam::Time::timeName()
  // cf. src/OpenFOAM/db/Time/Time.C
  vtksys_ios::ostringstream parser;
#ifdef _MSC_VER
  bool correctExponent = true;
#endif
  if(timeFormat == "general")
    {
    // "do not use std:: or vtkstd:: when using anything declared in
    // iostream," according to VTK coding standards, but following
    // the instruction causes errors...
    parser.setf(vtkstd::ios_base::fmtflags(0), vtkstd::ios_base::floatfield);
    }
  else if(timeFormat == "fixed")
    {
    parser.setf(vtkstd::ios_base::fmtflags(vtkstd::ios_base::fixed),
      vtkstd::ios_base::floatfield);
#ifdef _MSC_VER
    correctExponent = false;
#endif
    }
  else if(timeFormat == "scientific")
    {
    parser.setf(vtkstd::ios_base::fmtflags(vtkstd::ios_base::scientific),
      vtkstd::ios_base::floatfield);
    }
  else
    {
    vtkWarningMacro("Warning: unsupported time format. Assuming general.");
    parser.setf(vtkstd::ios_base::fmtflags(0), vtkstd::ios_base::floatfield);
    }
  parser.precision(timePrecision);

  for(int i = 0; i < tempNumTimeSteps; i++)
    {
    parser.str("");
    const double tempStep = i * timeStepIncrement + startTime;
    parser << tempStep; // stringstream doesn't require ends
#ifdef _MSC_VER
    // workaround for format difference in MSVC++:
    // remove an extra 0 from exponent
    if(correctExponent)
      {
      vtkStdString tempStr(parser.str());
      vtkStdString::size_type pos = tempStr.find('e');
      if(pos != vtkStdString::npos && tempStr.length() >= pos + 3
        && tempStr[pos + 2] == '0')
        {
        tempStr.erase(pos + 2, 1);
        parser.str(tempStr);
        }
      }
#endif
    // Add the time steps that actually exist to steps
    // allows the run to be stopped short of controlDict spec
    // allows for removal of timesteps
    if(test->Open((this->CasePath + parser.str()).c_str()))
      {
      this->TimeValues->InsertNextValue(tempStep);
      this->TimeNames->InsertNextValue(parser.str());
      }
    // necessary for reading the case/0 directory whatever the timeFormat is
    // based on Foam::Time::operator++() cf. src/OpenFOAM/db/Time/Time.C
    else if((fabs(tempStep) < 1.0e-14L) // 10*SMALL
      && test->Open((this->CasePath + vtkStdString("0")).c_str()))
      {
      this->TimeValues->InsertNextValue(tempStep);
      this->TimeNames->InsertNextValue(vtkStdString("0"));
      }
    }
  test->Delete();
  this->TimeValues->Squeeze();
  this->TimeNames->Squeeze();

  if(this->TimeValues->GetNumberOfTuples() == 0)
    {
    // set the number of timesteps to 1 if the constant subdirectory exists
    test = vtkDirectory::New();
    if(test->Open((this->CasePath + "constant").c_str()))
      {
      parser.str("");
      parser << startTime;
      this->TimeValues->InsertNextValue(startTime);
      this->TimeValues->Squeeze();
      this->TimeNames->InsertNextValue(parser.str());
      this->TimeNames->Squeeze();
      }
    test->Delete();
    }
  return true;
}

//-----------------------------------------------------------------------------
// list time directories by searching all valid time instances in a
// case directory
bool vtkOpenFOAMReaderPrivate::ListTimeDirectoriesByInstances()
{
  // open the case directory
  vtkDirectory* test = vtkDirectory::New();
  if(!test->Open(this->CasePath.c_str()))
    {
    test->Delete();
    vtkErrorMacro(<< "Can't open directory " << this->CasePath.c_str());
    return false;
    }

  // search all the directories in the case directory and detect
  // directories with names convertible to numbers
  this->TimeValues->Initialize();
  this->TimeNames->Initialize();
  const int nFiles = test->GetNumberOfFiles();
  for(int i = 0; i < nFiles; i++)
    {
    const vtkStdString dir = test->GetFile(i);
    if(test->FileIsDirectory(dir.c_str()))
      {
      // check if the name is convertible to a number
      bool isTimeDir = true;
      for(size_t j = 0; j < dir.length(); j++)
        {
        const char c = dir[j];
        if(!isdigit(c) && c != '+' && c != '-' && c != '.' && c != 'e'
          && c != 'E')
          {
          isTimeDir = false;
          break;
          }
        }
      if(!isTimeDir)
        {
        continue;
        }

      // convert to a number
      char *endptr;
      double timeValue = strtod(dir.c_str(), &endptr);
      // check if the value really was converted to a number
      if(timeValue == 0.0 && endptr == dir.c_str())
        {
        continue;
        }

      // add to the instance list
      this->TimeValues->InsertNextValue(timeValue);
      this->TimeNames->InsertNextValue(dir);
      }
    }
  test->Delete();
  this->TimeValues->Squeeze();
  this->TimeNames->Squeeze();

  // sort the detected time directories and set as timesteps
  if(this->TimeValues->GetNumberOfTuples() > 1)
    {
    vtkSortDataArray::Sort(this->TimeValues, this->TimeNames);
    }
  else if(this->TimeValues->GetNumberOfTuples() == 0)
    {
    // set the number of timesteps to 1 if the constant subdirectory exists
    test = vtkDirectory::New();
    if(test->Open((this->CasePath + "constant").c_str()))
      {
      this->TimeValues->InsertNextValue(0.0);
      this->TimeValues->Squeeze();
      this->TimeNames->InsertNextValue("0");
      this->TimeNames->Squeeze();
      }
    test->Delete();
    }

  return true;
}

//-----------------------------------------------------------------------------
// gather the necessary information to create a path to the data
bool vtkOpenFOAMReaderPrivate::MakeInformationVector(
  const vtkStdString &casePath, const vtkStdString &controlDictPath,
  const vtkStdString &procName, vtkOpenFOAMReader *parent)
{
  this->CasePath = casePath;
  this->ProcessorName = procName;
  this->Parent = parent;
  
  // list timesteps (skip parsing controlDict entirely if
  // ListTimeStepsByControlDict is set to 0)
  bool ret;
  if(this->Parent->GetListTimeStepsByControlDict())
    {
    vtkFoamIOobject io(this->CasePath);

    // open and check if controlDict is readable
    if(!io.open(controlDictPath))
      {
      vtkErrorMacro(<<"Error opening " << io.fileName().c_str() << ": "
        << io.error().c_str());
      return false;
      }
    vtkFoamDict dict;
    if(!dict.read(io))
      {
      vtkErrorMacro(<<"Error reading line " << io.lineNumber()
        << " of " << io.fileName().c_str() << ": " << io.error().c_str());
      return false;
      }
    if(dict.type() != vtkFoamToken::DICTIONARY)
      {
      vtkErrorMacro(<<"The file type of " << io.fileName().c_str()
        << " is not a dictionary");
      return false;
      }

    vtkFoamEntry& writeControlEntry = dict.lookup("writeControl");
    if(!writeControlEntry.found())
      {
      vtkErrorMacro(<< "writeControl entry not found in "
        << io.fileName().c_str());
      return false;
      }
    const vtkStdString writeControl(writeControlEntry.toString());

    // empty if not found
    const vtkStdString adjustTimeStep(dict.lookup("adjustTimeStep").toString());

    // list time directories according to controlDict if (adjustTimeStep
    // writeControl) == (off, timeStep) or (on, adjustableRunTime); list
    // by time instances in the case directory otherwise (different behavior
    // from paraFoam)
    // valid switching words cf. src/OpenFOAM/db/Switch/Switch.C
    if((((adjustTimeStep == "off" || adjustTimeStep == "no"
      || adjustTimeStep == "n" || adjustTimeStep == "false"
      || adjustTimeStep == "") && writeControl == "timeStep")
      || ((adjustTimeStep == "on" || adjustTimeStep == "yes"
      || adjustTimeStep == "y" || adjustTimeStep == "true")
      && writeControl == "adjustableRunTime")))
      {
      ret = this->ListTimeDirectoriesByControlDict(&dict);
      }
    else
      {
      ret = this->ListTimeDirectoriesByInstances();
      }
    }
  else
    {
    ret = this->ListTimeDirectoriesByInstances();
    }

  if(!ret)
    {
    return ret;
    }

  if(this->TimeStep >= this->TimeValues->GetNumberOfTuples())
    {
    this->SetTimeStep(this->TimeValues->GetNumberOfTuples() - 1);
    }

  this->PopulatePolyMeshDirArrays();
  return ret;
}

//-----------------------------------------------------------------------------
void vtkOpenFOAMReaderPrivate::AppendMeshDirToArray(vtkStringArray* polyMeshDir,
  const vtkStdString &path, const int timeI)
{
  vtkFoamIOobject io(this->CasePath);

  if(io.open(path) || io.open(path + ".gz"))
    {
    io.close();
    // set points/faces location to current timesteps value
    polyMeshDir->SetValue(timeI, this->TimeNames->GetValue(timeI));
    }
  else
    {
    if(timeI != 0)
      {
      // set points/faces location to previous timesteps value
      polyMeshDir->SetValue(timeI, polyMeshDir->GetValue(timeI - 1));
      }
    else
      {
      // set points/faces to constant
      polyMeshDir->SetValue(timeI, "constant");
      }
    }
}

//-----------------------------------------------------------------------------
// create a Lookup Table containing the location of the points
// and faces files for each time steps mesh
void vtkOpenFOAMReaderPrivate::PopulatePolyMeshDirArrays()
{
  // intialize size to number of timesteps
  const int nSteps = this->TimeValues->GetNumberOfTuples();
  this->PolyMeshPointsDir->SetNumberOfValues(nSteps);
  this->PolyMeshFacesDir->SetNumberOfValues(nSteps);

  // loop through each timestep
  for(int i = 0; i < nSteps; i++)
    {
    // create the path to the timestep
    vtkStdString polyMeshPath = this->TimeRegionPath(i) + "/polyMesh/";
    AppendMeshDirToArray(this->PolyMeshPointsDir, polyMeshPath + "points", i);
    AppendMeshDirToArray(this->PolyMeshFacesDir, polyMeshPath + "faces", i);
    }
  return;
}

//-----------------------------------------------------------------------------
// read the points file into a vtkFloatArray
vtkFloatArray* vtkOpenFOAMReaderPrivate::ReadPointsFile()
{
  // path to points file
  const vtkStdString pointPath
    = this->CurrentTimeRegionMeshPath(this->PolyMeshPointsDir) + "points";

  vtkFoamIOobject io(this->CasePath);
  if(!(io.open(pointPath) || io.open(pointPath + ".gz")))
    {
    vtkErrorMacro(<<"Error opening " << io.fileName().c_str() << ": "
      << io.error().c_str());
    return NULL;
    }

  vtkFoamEntryValue dict(NULL);
  try
    {
    dict.readNonuniformList<vtkFoamToken::VECTORLIST,
      vtkFoamEntryValue::vectorListTraits<vtkFloatArray, float, 3, false> >(io);
    }
  catch(vtkFoamError& e)
    {
    vtkErrorMacro(<<"Error reading line " << io.lineNumber()
      << " of " << io.fileName().c_str() << ": " << e.c_str());
    return NULL;
    }

  vtkFloatArray *pointArray = static_cast<vtkFloatArray *>(dict.ptr());

  // set the number of points
  this->NumPoints = pointArray->GetNumberOfTuples();

  return pointArray;
}

//-----------------------------------------------------------------------------
// read the faces into a intVectorVector
vtkOpenFOAMReaderPrivate::intVectorVector *
  vtkOpenFOAMReaderPrivate::ReadFacesFile(const vtkStdString &facePathIn)
{
  const vtkStdString facePath(facePathIn + "faces");

  vtkFoamIOobject io(this->CasePath);
  if(!(io.open(facePath) || io.open(facePath + ".gz")))
    {
    vtkErrorMacro(<<"Error opening " << io.fileName().c_str() << ": "
      << io.error().c_str());
    return NULL;
    }

  vtkFoamEntryValue dict(NULL);
  try
    {
    dict.readLabelListList(io);
    }
  catch(vtkFoamError& e)
    {
    vtkErrorMacro(<<"Error reading line " << io.lineNumber()
      << " of " << io.fileName().c_str() << ": " << e.c_str());
    return NULL;
    }
  return static_cast<intVectorVector *>(dict.ptr());
}

//-----------------------------------------------------------------------------
// read the owner and neighbor file and create cellFaces
vtkOpenFOAMReaderPrivate::intVectorVector *
  vtkOpenFOAMReaderPrivate::ReadOwnerNeighborFiles(
  const vtkStdString &ownerNeighborPath, intVectorVector *facePoints)
{
  vtkFoamIOobject io(this->CasePath);
  vtkStdString ownerPath(ownerNeighborPath + "owner");
  if(io.open(ownerPath) || io.open(ownerPath + ".gz"))
    {
    vtkFoamEntryValue ownerDict(NULL);
    try
      {
      ownerDict.readNonuniformList<vtkFoamToken::LABELLIST,
        vtkFoamEntryValue::listTraits<vtkIntArray, int> >(io);
      }
    catch(vtkFoamError& e)
      {
      vtkErrorMacro(<<"Error reading line " << io.lineNumber()
      << " of " << io.fileName().c_str() << ": " << e.c_str());
      return NULL;
      }

    io.close();

    const vtkStdString neighborPath(ownerNeighborPath + "neighbour");
    if(!(io.open(neighborPath) || io.open(neighborPath + ".gz")))
      {
      vtkErrorMacro(<<"Error opening " << io.fileName().c_str() << ": "
      << io.error().c_str());
      return NULL;
      }

    vtkFoamEntryValue neighborDict(NULL);
    try
      {
      neighborDict.readNonuniformList<vtkFoamToken::LABELLIST,
        vtkFoamEntryValue::listTraits<vtkIntArray, int> >(io);
      }
    catch(vtkFoamError& e)
      {
      vtkErrorMacro(<<"Error reading line " << io.lineNumber()
      << " of " << io.fileName().c_str() << ": " << e.c_str());
      return NULL;
      }

    this->FaceOwner = static_cast<vtkIntArray *>(ownerDict.ptr());
    vtkIntArray &faceOwner = *this->FaceOwner;
    vtkIntArray &faceNeighbor = neighborDict.labelList();

    const int nFaces = faceOwner.GetNumberOfTuples();
    const int nNeiFaces = faceNeighbor.GetNumberOfTuples();

    if(nFaces < nNeiFaces)
      {
      vtkErrorMacro(<<"Numbers of owner faces " << nFaces
      << " must be equal or larger than number of neighbor faces "
      << nNeiFaces);
      return NULL;
      }

    if(nFaces != facePoints->nElements())
      {
      vtkWarningMacro(<<"Numbers of faces in faces "
        << facePoints->nElements() << " and owner " << nFaces
        << " does not match");
      return NULL;
      }

    // add the face numbers to the correct cell cf. Terry's code and
    // src/OpenFOAM/meshes/primitiveMesh/primitiveMeshCells.C
    // find the number of cells
    int nCells = -1;
    for(int faceI = 0; faceI < nNeiFaces; faceI++)
      {
      const int ownerCell = faceOwner.GetValue(faceI);
      if(nCells < ownerCell) // max(nCells, faceOwner[i])
        {
        nCells = ownerCell;
        }
      // we do need to take neighbor faces into account since all the
      // surrounding faces of a cell can be neighbors for a valid mesh
      const int neighborCell = faceNeighbor.GetValue(faceI);
      if(nCells < neighborCell) // max(nCells, faceNeighbor[i])
        {
        nCells = neighborCell;
        }
      }
    for(int faceI = nNeiFaces; faceI < nFaces; faceI++)
      {
      const int ownerCell = faceOwner.GetValue(faceI);
      if(nCells < ownerCell) // max(nCells, faceOwner[i])
        {
        nCells = ownerCell;
        }
      }
    nCells++;

    if(nCells == 0)
      {
      vtkWarningMacro(<<"The mesh contains no cells");
      }

    // set the number of cells
    this->NumCells = nCells;

    // create cellFaces with the length of the body undetermined
    intVectorVector *cells = new intVectorVector(nCells, 1);

    // count number of faces for each cell
    int *cfiPtr = cells->indices()->GetPointer(0);
    for(int cellI = 0; cellI <= nCells; cellI++)
      {
      cfiPtr[cellI] = 0;
      }
    int nTotalCellFaces = 0;
    cfiPtr++; // offset +1
    for(int faceI = 0; faceI < nNeiFaces; faceI++)
      {
      const int ownerCell = faceOwner.GetValue(faceI);
      // simpleFoam/pitzDaily3Blocks has faces with owner cell number -1
      if(ownerCell >= 0)
        {
        cfiPtr[ownerCell]++;
        nTotalCellFaces++;
        }
      const int neighborCell=faceNeighbor.GetValue(faceI);
      if(neighborCell >= 0)
        {
        cfiPtr[neighborCell]++;
        nTotalCellFaces++;
        }
      }
    for(int faceI = nNeiFaces; faceI < nFaces; faceI++)
      {
      const int ownerCell = faceOwner.GetValue(faceI);
      if(ownerCell >= 0)
        {
        cfiPtr[ownerCell]++;
        nTotalCellFaces++;
        }
      }
    cfiPtr--; // revert offset +1

    // allocate cellFaces. To reduce the numbers of new/delete operations we
    // allocate memory space for all faces linearly
    cells->resizeBody(nTotalCellFaces);

    // accumulate the number of cellFaces to create cellFaces indices
    // and copy them to a temporary array
    vtkIntArray *tmpFaceIndices = vtkIntArray::New();
    tmpFaceIndices->SetNumberOfValues(nCells + 1);
    int *tfiPtr = tmpFaceIndices->GetPointer(0);
    tfiPtr[0] = 0;
    for(int cellI = 1; cellI <= nCells; cellI++)
      {
      tfiPtr[cellI] = (cfiPtr[cellI] += cfiPtr[cellI - 1]);
      }

    // add face numbers to cell-faces list
    vtkIntArray *cellFacesList = cells->body();
    for(int faceI = 0; faceI < nNeiFaces; faceI++)
      {
      const int ownerCell = faceOwner.GetValue(faceI); // must be a signed int
      // simpleFoam/pitzDaily3Blocks has faces with owner cell number -1
      if(ownerCell >= 0)
        {
        cellFacesList->SetValue(tfiPtr[ownerCell]++, faceI);
        }
      const int neighborCell = faceNeighbor.GetValue(faceI);
      if(neighborCell >= 0)
        {
        cellFacesList->SetValue(tfiPtr[neighborCell]++, faceI);
        }
      }
    for(int faceI = nNeiFaces; faceI < nFaces; faceI++)
      {
      const int ownerCell = faceOwner.GetValue(faceI); // must be a signed int
      // simpleFoam/pitzDaily3Blocks has faces with owner cell number -1
      if(ownerCell >= 0)
        {
        cellFacesList->SetValue(tfiPtr[ownerCell]++, faceI);
        }
      }
    tmpFaceIndices->Delete();

    return cells;
    }
  else // if owner does not exist look for cells
    {
    vtkStdString cellsPath(ownerNeighborPath + "cells");
    if(!(io.open(cellsPath) || io.open(cellsPath + ".gz")))
      {
      vtkErrorMacro(<<"Error opening " << io.fileName().c_str() << ": "
        << io.error().c_str());
      return NULL;
      }
    vtkFoamEntryValue cellsDict(NULL);
    try
      {
      cellsDict.readLabelListList(io);
      }
    catch(vtkFoamError& e)
      {
      vtkErrorMacro(<<"Error reading line " << io.lineNumber()
        << " of " << io.fileName().c_str() << ": " << e.c_str());
      return NULL;
      }

    intVectorVector *cells = static_cast<intVectorVector *>(cellsDict.ptr());
    this->NumCells = cells->nElements();
    const int nFaces = facePoints->nElements();

    // create face owner list
    this->FaceOwner = vtkIntArray::New();
    this->FaceOwner->SetNumberOfTuples(nFaces);
    for(int faceI = 0; faceI < nFaces; faceI++)
      {
      this->FaceOwner->SetValue(faceI, -1);
      }
    for(int cellI = 0; cellI < this->NumCells; cellI++)
      {
      const int nCellFaces = cells->size(cellI);
      const int *cellFaces = cells->operator[](cellI);
      for(int faceI = 0; faceI < nCellFaces; faceI++)
        {
        const int f = cellFaces[faceI];
        if(f < 0 || f >= nFaces) // make sure the face number is valid
          {
          vtkErrorMacro("Face number " << f << " in cell " << cellI
          << " exceeds the number of faces " << nFaces);
          this->FaceOwner->Delete();
          this->FaceOwner = NULL;
          delete cells;
          return NULL;
          }
        const int owner = this->FaceOwner->GetValue(f);
        if(owner == -1 || owner > cellI)
          {
          this->FaceOwner->SetValue(f, cellI);
          }
        }
      }
    // check for unused faces
    for(int faceI = 0; faceI < nFaces; faceI++)
      {
      if(this->FaceOwner->GetValue(faceI) == -1)
        {
        vtkErrorMacro(<<"Face " << faceI << " is not used");
        this->FaceOwner->Delete();
        this->FaceOwner = NULL;
        delete cells;
        return NULL;
        }
      }
    return cells;
    }
}

//-----------------------------------------------------------------------------
bool vtkOpenFOAMReaderPrivate::CheckFacePoints(intVectorVector *facePoints)
{
  const int nFaces = facePoints->nElements();

  for(int faceI = 0; faceI < nFaces; faceI++)
    {
    const int nPoints = facePoints->size(faceI);
    const int *pointList = facePoints->operator[](faceI);
    for(int pointI = 0; pointI < nPoints; pointI++)
      {
      const int p = pointList[pointI];
      if(p < 0 || p >= this->NumPoints)
        {
        vtkErrorMacro(<< "The point number " << p << " at face number " << faceI
          << " is out of range for " << this->NumPoints << " points");
        return false;
        }
      }
    }
  return true;
}

//-----------------------------------------------------------------------------
// determine cell shape and insert the cell into the mesh
// hexahedron, prism, pyramid, tetrahedron and decompose polyhedron
void vtkOpenFOAMReaderPrivate::InsertCellsToGrid(
  vtkUnstructuredGrid* internalMesh, const intVectorVector *cellsFaces,
  const intVectorVector *facesPoints, vtkFloatArray *pointArray,
  vtkIdTypeArray *additionalCells, vtkIntArray *cellList)
{
  const int maxNPoints = 128; // assume max number of points per cell
  vtkIdList* cellPoints = vtkIdList::New();
  cellPoints->SetNumberOfIds(maxNPoints);
  const int nCells
    = (cellList == NULL ? this->NumCells : cellList->GetNumberOfTuples());
  int nAdditionalPoints = 0;

  // alias
  const intVectorVector& facePoints = *facesPoints;

  for(int cellI = 0; cellI < nCells ; cellI++)
    {
    int cellId;
    if(cellList == NULL)
      {
      cellId = cellI;
      }
    else
      {
      cellId = cellList->GetValue(cellI);
      if(cellId >= this->NumCells)
        {
        vtkWarningMacro(<<"cellLabels id " << cellId
          << " exceeds the number of cells " << nCells
          << ". Inserting an empty cell.");
        internalMesh->InsertNextCell(VTK_EMPTY_CELL, 0,
          cellPoints->GetPointer(0));
        continue;
        }
      }
    const int *cellFaces = cellsFaces->operator[](cellId);
    const int nCellFaces = cellsFaces->size(cellId);

    // determine type of the cell
    // cf. src/OpenFOAM/meshes/meshShapes/cellMatcher/{hex|prism|pyr|tet}-
    // Matcher.C
    int cellType = VTK_CONVEX_POINT_SET;
    if(nCellFaces == 6)
      {
      int j = 0;
      for(; j < nCellFaces; j++)
        {
        if(facePoints.size(cellFaces[j]) != 4)
          {
          break;
          }
        }
      if(j == nCellFaces)
        {
        cellType = VTK_HEXAHEDRON;
        }
      }
    else if(nCellFaces == 5)
      {
      int nTris = 0, nQuads = 0;
      for(int j = 0; j < nCellFaces; j++)
        {
        const int nPoints = facePoints.size(cellFaces[j]);
        if(nPoints == 3)
          {
          nTris++;
          }
        else if(nPoints == 4)
          {
          nQuads++;
          }
        else
          {
          break;
          }
        }
      if(nTris == 2 && nQuads == 3)
        {
        cellType = VTK_WEDGE;
        }
      else if(nTris == 4 && nQuads == 1)
        {
        cellType = VTK_PYRAMID;
        }
      }
    else if(nCellFaces == 4)
      {
      int j = 0;
      for(; j < nCellFaces; j++)
        {
        if(facePoints.size(cellFaces[j]) != 3)
          {
          break;
          }
        }
      if(j == nCellFaces)
        {
        cellType = VTK_TETRA;
        }
      }

    // not a Hex/Wedge/Pyramid/Tetra
    if(cellType == VTK_CONVEX_POINT_SET)
      {
      int nPoints = 0;
      for(int j = 0; j < nCellFaces; j++)
        {
        nPoints += facePoints.size(cellFaces[j]);
        }
      if(nPoints == 0)
        {
        cellType = VTK_EMPTY_CELL;
        }
      }

    // Cell shape constructor based on the one implementd by Terry
    // Jordan, with lots of improvements. Not as elegant as the one in
    // OpenFOAM but it's simple and works reasonably fast.

    // OFhex | vtkHexahedron
    if (cellType == VTK_HEXAHEDRON)
      {
      // get first face in correct order
      const int cellBaseFaceId = cellFaces[0];
      const int *face0Points = facePoints[cellBaseFaceId];

      if(this->FaceOwner->GetValue(cellBaseFaceId) == cellId)
        {
        // if it is an owner face flip the points
        for(int j = 0; j < 4; j++)
          {
          cellPoints->SetId(j, face0Points[3 - j]);
          }
        }
      else
        {
        // add base face to cell points
        for(int j = 0; j < 4; j++)
          {
          cellPoints->SetId(j, face0Points[j]);
          }
        }
      const int baseFacePoint0 = cellPoints->GetId(0);
      const int baseFacePoint2 = cellPoints->GetId(2);
      int cellOppositeFaceI = -1, pivotPoint = -1;
      int dupPoint = -1;
      for(int faceI = 1; faceI < 5; faceI++) // skip face 0 and 5
        {
        const int cellFaceI = cellFaces[faceI];
        const int *faceIPoints = facePoints[cellFaceI];
        int foundDup = -1, pointI = 0;
        for(; pointI < 4; pointI++) // each point
          {
          const int faceIPointI = faceIPoints[pointI];
          // matching two points in base face is enough to find a
          // duplicated point since neighboring faces share two
          // neighboring points (i. e. an edge)
          if(baseFacePoint0 == faceIPointI)
            {
            foundDup = 0;
            break;
            }
          else if(baseFacePoint2 == faceIPointI)
            {
            foundDup = 2;
            break;
            }
          }
        if(foundDup >= 0)
          {
          // find the pivot point if still haven't
          if(pivotPoint == -1)
            {
            dupPoint = foundDup;

            const int faceINextPoint = faceIPoints[(pointI + 1) % 4];

            // if the next point of the faceI-th face matches the
            // previous point of the base face use the previous point
            // of the faceI-th face as the pivot point; or use the
            // next point otherwise
            if(faceINextPoint == (this->FaceOwner->GetValue(cellFaceI)
              == cellId ? cellPoints->GetId(1 + foundDup)
              : cellPoints->GetId(3 - foundDup)))
              {
              pivotPoint = faceIPoints[(3 + pointI) % 4];
              }
            else
              {
              pivotPoint = faceINextPoint;
              }

            if(cellOppositeFaceI >= 0)
              {
              break;
              }
            }
          }
        else
          {
          // if no duplicated point found, faceI is the opposite face
          cellOppositeFaceI = cellFaceI;

          if(pivotPoint >= 0)
            {
            break;
            }
          }
        }

      // if the opposite face is not found until face 4, face 5 is
      // always the opposite face
      if(cellOppositeFaceI == -1)
        {
        cellOppositeFaceI = cellFaces[5];
        }

      // find the pivot point in opposite face
      const int *oppositeFacePoints = facePoints[cellOppositeFaceI];
      int pivotPointI = 0;
      for(; pivotPointI < 4; pivotPointI++)
        {
        if(oppositeFacePoints[pivotPointI] == pivotPoint)
          {
          break;
          }
        }

      // shift the pivot point if the point corresponds to point 2
      // of the base face
      if(dupPoint == 2)
	{
        pivotPointI = (pivotPointI + 2) % 4;
        }
      // copy the face-point list of the opposite face to cell-point list
      int basePointI = 4;
      if(this->FaceOwner->GetValue(cellOppositeFaceI) == cellId)
        {
        for(int pointI = pivotPointI; pointI < 4; pointI++)
          {
          cellPoints->SetId(basePointI++, oppositeFacePoints[pointI]);
          }
        for(int pointI = 0; pointI < pivotPointI; pointI++)
          {
          cellPoints->SetId(basePointI++, oppositeFacePoints[pointI]);
          }
        }
      else
        {
        for(int pointI = pivotPointI; pointI >= 0; pointI--)
          {
          cellPoints->SetId(basePointI++, oppositeFacePoints[pointI]);
          }
        for(int pointI = 3; pointI > pivotPointI; pointI--)
          {
          cellPoints->SetId(basePointI++, oppositeFacePoints[pointI]);
          }
        }

      // create the hex cell and insert it into the mesh
      internalMesh->InsertNextCell(cellType, 8, cellPoints->GetPointer(0));
      }

    // the cell construction is about the same as that of a hex, but
    // the point ordering have to be reversed!!
    else if (cellType == VTK_WEDGE)
      {
      // find the base face number
      int baseFaceId = 0;
      for(int j = 0; j < 5; j++)
        {
        if(facePoints.size(cellFaces[j]) == 3)
          {
          baseFaceId = j;
          break;
          }
        }

      // get first face in correct order
      const int cellBaseFaceId = cellFaces[baseFaceId];
      const int *face0Points = facePoints[cellBaseFaceId];

      if(this->FaceOwner->GetValue(cellBaseFaceId) == cellId)
        {
        for(int j = 0; j < 3; j++)
          {
          cellPoints->SetId(j, face0Points[j]);
          }
        }
      else
        {
        // if it is a neighbor face flip the points
        for(int j = 0; j < 3; j++)
          {
          // add base face to cell points
          cellPoints->SetId(j, face0Points[2 - j]);
          }
        }
      const int baseFacePoint0 = cellPoints->GetId(0);
      const int baseFacePoint2 = cellPoints->GetId(2);
      int cellOppositeFaceI = -1, pivotPoint = -1;
      bool dupPoint2 = false;
      for(int faceI = 0; faceI < 5; faceI++)
        {
        if(faceI == baseFaceId)
          {
          continue;
          }
        const int cellFaceI = cellFaces[faceI];
        if(facePoints.size(cellFaceI) == 3)
          {
          cellOppositeFaceI = cellFaceI;
          }
        // find the pivot point if still haven't
        else if(pivotPoint == -1)
          {
          const int *faceIPoints = facePoints[cellFaceI];
          bool found0Dup = false, found2Dup = false;
          int pointI = 0;
          for(; pointI < 4; pointI++) // each point
            {
            const int faceIPointI = faceIPoints[pointI];
            // matching two points in base face is enough to find a
            // duplicated point since neighboring faces share two
            // neighboring points (i. e. an edge)
            if(baseFacePoint0 == faceIPointI)
              {
              found0Dup = true;
              break;
              }
            else if(baseFacePoint2 == faceIPointI)
              {
              found2Dup = true;
              break;
              }
            }
          // the matching point must always be found so omit the check
          int baseFacePrevPoint, baseFaceNextPoint;
          if(found0Dup)
            {
            baseFacePrevPoint = cellPoints->GetId(2);
            baseFaceNextPoint = cellPoints->GetId(1);
            }
          else
            {
            baseFacePrevPoint = cellPoints->GetId(1);
            baseFaceNextPoint = cellPoints->GetId(0);
            dupPoint2 = true;
            }

          const int faceINextPoint = faceIPoints[(pointI + 1) % 4];
          const int faceIPrevPoint = faceIPoints[(3 + pointI) % 4];

          // if the next point of the faceI-th face matches the
          // previous point of the base face use the previous point of
          // the faceI-th face as the pivot point; or use the next
          // point otherwise
          if(faceINextPoint == (this->FaceOwner->GetValue(cellFaceI)
            == cellId ? baseFacePrevPoint : baseFaceNextPoint))
            {
            pivotPoint = faceIPrevPoint;
            }
          else
            {
            pivotPoint = faceINextPoint;
            }
          }

        // break when both of opposite face and pivot point are found
        if(cellOppositeFaceI >= 0 && pivotPoint >= 0)
          {
          break;
          }
        }

      // find the pivot point in opposite face
      const int *oppositeFacePoints = facePoints[cellOppositeFaceI];
      int pivotPointI = 0;
      for(; pivotPointI < 3; pivotPointI++)
        {
        if(oppositeFacePoints[pivotPointI] == pivotPoint)
          {
          break;
          }
        }

      if(this->FaceOwner->GetValue(cellOppositeFaceI) == cellId)
        {
        if(dupPoint2)
          {
          pivotPointI = (pivotPointI + 2) % 3;
          }
        int basePointI = 3;
        for(int pointI = pivotPointI; pointI >= 0; pointI--)
          {
          cellPoints->SetId(basePointI++, oppositeFacePoints[pointI]);
          }
        for(int pointI = 2; pointI > pivotPointI; pointI--)
          {
          cellPoints->SetId(basePointI++, oppositeFacePoints[pointI]);
          }
        }
      else
        {
        // shift the pivot point if the point corresponds to point 2
        // of the base face
        if(dupPoint2)
          {
          pivotPointI = (1 + pivotPointI) % 3;
          }
        // copy the face-point list of the opposite face to cell-point list
        int basePointI = 3;
        for(int pointI = pivotPointI; pointI < 3; pointI++)
          {
          cellPoints->SetId(basePointI++, oppositeFacePoints[pointI]);
          }
        for(int pointI = 0; pointI < pivotPointI; pointI++)
          {
          cellPoints->SetId(basePointI++, oppositeFacePoints[pointI]);
          }
        }

      // create the wedge cell and insert it into the mesh
      internalMesh->InsertNextCell(cellType, 6, cellPoints->GetPointer(0));
      }

    // OFpyramid | vtkPyramid || OFtet | vtkTetrahedron
    else if (cellType == VTK_PYRAMID || cellType == VTK_TETRA)
      {
      int baseFaceId = -1, nPoints;
      if(cellType == VTK_PYRAMID)
        {
        for(int j = 0; j < nCellFaces; j++)
          {
          if(facePoints.size(cellFaces[j]) == 4)
            {
            baseFaceId = j;
            break;
            }
          }
        nPoints = 5;
        }
      else // VTK_TETRA
        {
        baseFaceId = 0;
        nPoints = 4;
        }

      // add first face to cell points
      const int cellBaseFaceId = cellFaces[baseFaceId];
      const int *baseFacePoints = facePoints[cellBaseFaceId];
      const size_t nBaseFacePoints = facePoints.size(cellBaseFaceId);
      if(this->FaceOwner->GetValue(cellBaseFaceId) == cellId)
        {
        // if it is an owner face flip the points
        for(size_t j = 0; j < nBaseFacePoints; j++)
          {
          cellPoints->SetId(j, baseFacePoints[nBaseFacePoints - 1 - j]);
          }
        }
      else
        {
        for(size_t j = 0; j < nBaseFacePoints; j++)
          {
          cellPoints->SetId(j, baseFacePoints[j]);
          }
        }

      // compare an adjacent face (any non base face is ok) point 1 to
      // base face points
      const int adjacentFaceId = (baseFaceId == 0) ? 1 : baseFaceId - 1;
      const int cellAdjacentFaceId = cellFaces[adjacentFaceId];
      const int *adjacentFacePoints = facePoints[cellAdjacentFaceId];
      const int adjacentFacePoint1 = adjacentFacePoints[1];
      bool foundDup = false;
      for(size_t j = 0; j < nBaseFacePoints; j++)
        {
        // if point 1 of the adjacent face matches point j of the base face...
        if(cellPoints->GetId(j) == adjacentFacePoint1)
          {
          // if point 2 of the adjacent face matches the previous point
          // of the base face use point 0 of the adjacent face as the
          // pivot point; use point 2 otherwise
          cellPoints->SetId(nBaseFacePoints, (adjacentFacePoints[2]
            == cellPoints->GetId((this->FaceOwner->GetValue(cellAdjacentFaceId)
            == cellId ? (j + 1) : (nBaseFacePoints + j - 1)) % nBaseFacePoints))
            ? adjacentFacePoints[0] : adjacentFacePoints[2]);
          foundDup = true;
          break;
          }
        }
      // if point 1 of the adjacent face does not match any points of
      // the base face, it's the pivot point
      if(!foundDup)
        {
        cellPoints->SetId(nBaseFacePoints, adjacentFacePoint1);
        }

      // create the tetra cell and insert it into the mesh
      internalMesh->InsertNextCell(cellType, nPoints,
        cellPoints->GetPointer(0));
      }

    // erronous cells
    else if(cellType == VTK_EMPTY_CELL)
      {
      vtkWarningMacro("Warning: No points in cellId " << cellId);
      internalMesh->InsertNextCell(VTK_EMPTY_CELL, 0,
        cellPoints->GetPointer(0));
      }

    // OFpolyhedron || vtkConvexPointSet
    else
      {
      if(additionalCells != NULL) // decompose into tets and pyramids
        {
        // calculate cell centroid and insert it to point list
        this->AdditionalCellPoints->push_back(vtkIntArray::New());
        vtkIntArray *polyCellPoints = this->AdditionalCellPoints->back();
        float centroid[3];
        centroid[0] = centroid[1] = centroid[2] = 0.0F;
        for(int j = 0; j < nCellFaces; j++)
          {
          // remove duplicate points from faces
          const int cellFacesJ = cellFaces[j];
          const int *faceJPoints = facePoints[cellFacesJ];
          const size_t nFaceJPoints = facePoints.size(cellFacesJ);
          for(size_t k = 0; k < nFaceJPoints; k++)
            {
            const int faceJPointK = faceJPoints[k];
            bool foundDup = false;
            for(size_t l = 0; l < polyCellPoints->GetDataSize(); l++)
              {
              if(polyCellPoints->GetValue(l) == faceJPointK)
                {
                foundDup = true;
                break; // look no more
                }
              }
            if(!foundDup)
              {
              polyCellPoints->InsertNextValue(faceJPointK);
              float *pointK = pointArray->GetPointer(3 * faceJPointK);
              centroid[0] += pointK[0];
              centroid[1] += pointK[1];
              centroid[2] += pointK[2];
              }
            }
          }
        polyCellPoints->Squeeze();
        float weight = 1.0F / static_cast<float>(polyCellPoints->GetDataSize());
        centroid[0] *= weight;
        centroid[1] *= weight;
        centroid[2] *= weight;
        pointArray->InsertNextTuple(centroid);

        // polyhedron decomposition.
        // a tweaked algorithm based on applications/utilities/postProcessing/
        // graphics/PVFoamReader/vtkFoam/vtkFoamAddInternalMesh.C
        bool insertDecomposedCell = true;
        for(int j = 0; j < nCellFaces; j++)
          {
          const int cellFacesJ = cellFaces[j];
          const int *faceJPoints = facePoints[cellFacesJ];
          const int nFaceJPoints = facePoints.size(cellFacesJ);
          const int flipNeighbor
            = (this->FaceOwner->GetValue(cellFacesJ) == cellId ? -1 : 1);
          const int nTris = nFaceJPoints % 2;

          int vertI = 2;

          // shift the start and end of the vertex loop if the
          // triangle of a decomposed face is going to be flat. Far
          // from perfect but better than nothing to avoid flat cells
          // which stops time integration of Stream Tracer especially
          // for split-hex unstructured meshes created by
          // e. g. autoRefineMesh
          if(nFaceJPoints >= 5 && nTris)
            {
            float *point0, *point1, *point2;
            point0 = pointArray->GetPointer(3 * faceJPoints[nFaceJPoints - 1]);
            point1 = pointArray->GetPointer(3 * faceJPoints[0]);
            point2 = pointArray->GetPointer(3 * faceJPoints[nFaceJPoints - 2]);
            float vsizeSqr1 = 0.0F, vsizeSqr2 = 0.0F, dotProduct = 0.0F;
            for(int i = 0; i < 3; i++)
              {
              const float v1 = point1[i] - point0[i],
                v2 = point2[i] - point0[i];
              vsizeSqr1 += v1 * v1;
              vsizeSqr2 += v2 * v2;
              dotProduct += v1 * v2;
              }
            // compare in squared representation to avoid using sqrt()
            if(dotProduct * fabsf(dotProduct) / (vsizeSqr1 * vsizeSqr2)
              < -1.0F + 1.0e-3F)
              {
              vertI = 1;
              }
            }

          cellPoints->SetId(0, faceJPoints[vertI == 2 ? 0 : nFaceJPoints - 1]);
          cellPoints->SetId(4, this->NumPoints + nAdditionalPoints);

          // decompose a face into quads in order (flipping the
          // decomposed face if owner)
          const int nQuadVerts = nFaceJPoints - 1 - nTris;
          for(; vertI < nQuadVerts; vertI += 2)
            {
            cellPoints->SetId(1, faceJPoints[vertI - flipNeighbor]);
            cellPoints->SetId(2, faceJPoints[vertI]);
            cellPoints->SetId(3, faceJPoints[vertI + flipNeighbor]);

            // if the decomposed cell is the first one insert it to
            // the original position; or append to the decomposed cell
            // list otherwise
            if(insertDecomposedCell)
              {
              internalMesh->InsertNextCell(VTK_PYRAMID, 5,
                cellPoints->GetPointer(0));
              insertDecomposedCell = false;
              }
            else
              {
              this->AdditionalCellIds->InsertNextValue(cellId);
              additionalCells->InsertNextTupleValue(cellPoints->GetPointer(0));
              }
            }

          // if the number of vertices is odd there's a triangle
          if(nTris)
            {
            if(flipNeighbor == -1)
              {
              cellPoints->SetId(1, faceJPoints[vertI]);
              cellPoints->SetId(2, faceJPoints[vertI - 1]);
              }
            else
              {
              cellPoints->SetId(1, faceJPoints[vertI - 1]);
              cellPoints->SetId(2, faceJPoints[vertI]);
              }
            cellPoints->SetId(3, this->NumPoints + nAdditionalPoints);

            if(insertDecomposedCell)
              {
              internalMesh->InsertNextCell(VTK_TETRA, 4,
                cellPoints->GetPointer(0));
              insertDecomposedCell = false;
              }
            else
              {
              // set the 5th vertex number to -1 to distinguish a tetra cell
              cellPoints->SetId(4, -1);
              this->AdditionalCellIds->InsertNextValue(cellId);
              additionalCells->InsertNextTupleValue(cellPoints->GetPointer(0));
              }
            }
          }
        nAdditionalPoints++;
        }
      else // don't decompose; use VTK_CONVEX_PONIT_SET
        {
        // get first face
        const int cellFaces0 = cellFaces[0];
        const int *baseFacePoints = facePoints[cellFaces0];
        const int nBaseFacePoints = facePoints.size(cellFaces0);
        int nPoints = nBaseFacePoints;
        if(nPoints > maxNPoints)
          {
          vtkErrorMacro(<< "Too large polyhedron at cellId = " << cellId);
          return;
          }
        if(this->FaceOwner->GetValue(cellFaces0) == cellId)
          {
          // if it is an owner face flip the points
          // not sure if flipping is necessary but do it anyway
          for(int j = 0; j < nBaseFacePoints; j++)
            {
            cellPoints->SetId(j, baseFacePoints[nBaseFacePoints - 1 - j]);
            }
          }
        else
          {
	  // add first face to cell points
          for(int j = 0; j < nBaseFacePoints; j++)
            {
            cellPoints->SetId(j, baseFacePoints[j]);
            }
          }

        // loop through faces and create a list of all points
        // j = 1 skip baseFace
        for(int j = 1; j < nCellFaces; j++)
          {
          // remove duplicate points from faces
          const int cellFacesJ = cellFaces[j];
          const int *faceJPoints = facePoints[cellFacesJ];
          const size_t nFaceJPoints = facePoints.size(cellFacesJ);
          for(size_t k = 0; k < nFaceJPoints; k++)
            {
            const int faceJPointK = faceJPoints[k];
            bool foundDup = false;
            for(int l = 0; l < nPoints; l++)
              {
              if(cellPoints->GetId(l) == faceJPointK)
                {
                foundDup = true;
                break; // look no more
                }
              }
            if(!foundDup)
              {
              if(nPoints >= maxNPoints)
                {
                vtkErrorMacro(<< "Too large polyhedron at cellId = " << cellId);
                return;
                }
              cellPoints->SetId(nPoints++, faceJPointK);
              }
            }
          }

        // create the poly cell and insert it into the mesh
        internalMesh->InsertNextCell(VTK_CONVEX_POINT_SET, nPoints,
          cellPoints->GetPointer(0));
        }
      }
    }
  cellPoints->Delete();
}

//-----------------------------------------------------------------------------
void vtkOpenFOAMReaderPrivate::SetBlockName(vtkMultiBlockDataSet *blocks,
  unsigned int blockI, const char *name)
{
  blocks->GetMetaData(blockI)->Set(vtkCompositeDataSet::NAME(), name);
}

//-----------------------------------------------------------------------------
// derive cell types and create the internal mesh
vtkUnstructuredGrid *vtkOpenFOAMReaderPrivate::MakeInternalMesh(
  const intVectorVector *cellsFaces, const intVectorVector *facesPoints,
  vtkFloatArray *pointArray)
{
  // Create Mesh
  vtkUnstructuredGrid* internalMesh = vtkUnstructuredGrid::New();
  internalMesh->Allocate(this->NumCells);

  if(this->Parent->GetDecomposePolyhedra())
    {
    // for polyhedral decomposition
    this->AdditionalCellIds = vtkIntArray::New();
    this->AdditionalCellPoints = new intArrayVector;

    vtkIdTypeArray *additionalCells = vtkIdTypeArray::New();
    additionalCells->SetNumberOfComponents(5); // accommodates tetra or pyramid

    this->InsertCellsToGrid(internalMesh, cellsFaces, facesPoints, pointArray,
      additionalCells, NULL);

    // for polyhedral decomposition
    pointArray->Squeeze();
    this->AdditionalCellIds->Squeeze();
    additionalCells->Squeeze();

    // insert decomposed cells into mesh
    const int nComponents = additionalCells->GetNumberOfComponents();
    const int nAdditionalCells = additionalCells->GetNumberOfTuples();
    for(int i = 0; i < nAdditionalCells; i++)
      {
      if(additionalCells->GetComponent(i, 4) == -1)
        {
        internalMesh->InsertNextCell(VTK_TETRA, 4,
          additionalCells->GetPointer(i * nComponents));
        }
      else
        {
        internalMesh->InsertNextCell(VTK_PYRAMID, 5,
          additionalCells->GetPointer(i * nComponents));
        }
      }
    internalMesh->Squeeze();
    additionalCells->Delete();
    }
  else
    {
    this->InsertCellsToGrid(internalMesh, cellsFaces, facesPoints, pointArray,
    NULL, NULL);
    }

  // set the internal mesh points
  vtkPoints *points = vtkPoints::New();
  points->SetData(pointArray);
  internalMesh->SetPoints(points);
  points->Delete();

  return internalMesh;
}

//-----------------------------------------------------------------------------
// insert faces to grid
void vtkOpenFOAMReaderPrivate::InsertFacesToGrid(
  vtkUnstructuredGrid *boundaryMesh, const intVectorVector *facesPoints,
  int startFace, int endFace, vtkIntArray *boundaryPointMap,
  vtkIdList *facePointsVtkId, vtkIntArray *labels, bool isLookupValue)
{
  vtkUnstructuredGrid &bm = *boundaryMesh;

  for(int j = startFace; j < endFace; j++)
    {
    int faceId;
    if(labels == NULL)
      {
      faceId = j;
      }
    else
      {
      faceId = labels->GetValue(j);
      if(faceId >= this->FaceOwner->GetNumberOfTuples())
        {
        vtkWarningMacro(<<"faceLabels id " << faceId
          << " exceeds the number of faces "
          << this->FaceOwner->GetNumberOfTuples());
        bm.InsertNextCell(VTK_EMPTY_CELL, 0, facePointsVtkId->GetPointer(0));
        continue;
        }
      }
    const int *facePoints = facesPoints->operator[](faceId);
    size_t nFacePoints = facesPoints->size(faceId);

    if(isLookupValue)
      {
      for(size_t k = 0; k < nFacePoints; k++)
        {
        facePointsVtkId->SetId(k, boundaryPointMap->LookupValue(facePoints[k]));
        }
      }
    else
      {
      if(boundaryPointMap)
        {
        for(size_t k = 0; k < nFacePoints; k++)
          {
          facePointsVtkId->SetId(k, boundaryPointMap->GetValue(facePoints[k]));
          }
        }
      else
        {
        for(size_t k = 0; k < nFacePoints; k++)
          {
          facePointsVtkId->SetId(k, facePoints[k]);
          }
        }
      }

    // triangle
    if(nFacePoints == 3)
      {
      bm.InsertNextCell(VTK_TRIANGLE, 3, facePointsVtkId->GetPointer(0));
      }
    // quad
    else if(nFacePoints == 4)
      {
      bm.InsertNextCell(VTK_QUAD, 4, facePointsVtkId->GetPointer(0));
      }
    // polygon
    else
      {
      bm.InsertNextCell(VTK_POLYGON, nFacePoints,
        facePointsVtkId->GetPointer(0));
      }
    }
}

//-----------------------------------------------------------------------------
// returns requested boundary meshes
vtkMultiBlockDataSet *vtkOpenFOAMReaderPrivate::MakeBoundaryMesh(
  const intVectorVector *facesPoints, vtkFloatArray* pointArray)
{
  const int nBoundaries = this->BoundaryDict.size();

  // do a consistency check of BoundaryDict
  int previousEndFace = -1;
  for(int boundaryI = 0; boundaryI < nBoundaries; boundaryI++)
    {
    const vtkFoamBoundaryEntry &beI = this->BoundaryDict[boundaryI];
    const int startFace = beI.startFace;
    const int nFaces = beI.nFaces;
    if(startFace < 0 || nFaces < 0)
      {
      vtkErrorMacro(<<"Neither of startFace " << startFace << " nor nFaces "
        << nFaces << " can be nagative for patch " << beI.boundaryName.c_str());
      return NULL;
      }
    if(previousEndFace >= 0 && previousEndFace != startFace)
      {
      vtkErrorMacro(<<"The end face number " << previousEndFace - 1
        << " of patch "
        << this->BoundaryDict[boundaryI - 1].boundaryName.c_str()
        << " is not consistent with the start face number " << startFace
        << " of patch " << beI.boundaryName.c_str());
      return NULL;
      }
    previousEndFace = startFace + nFaces;
    }
  if(previousEndFace != facesPoints->nElements())
    {
    vtkErrorMacro(<<"The end face number " << previousEndFace - 1
      << " of the last patch "
      << this->BoundaryDict[nBoundaries - 1].boundaryName.c_str()
      << " is not consistent with the number of faces "
      << facesPoints->nElements());
    return NULL;
    }

  vtkMultiBlockDataSet *boundaryMesh = vtkMultiBlockDataSet::New();

  if(this->Parent->GetCreateCellToPoint())
    {
    this->AllBoundaries = vtkUnstructuredGrid::New();
    this->AllBoundaries->Allocate(facesPoints->nElements()
      - this->BoundaryDict[0].startFace);
    }
  this->BoundaryPointMap = new intArrayVector;

  vtkIntArray *nBoundaryPointsList = vtkIntArray::New();
  nBoundaryPointsList->SetNumberOfValues(nBoundaries);

  // count the max number of points per face and the number of points
  // (with duplicates) in mesh
  int maxNFacePoints = 0;
  for(int boundaryI = 0; boundaryI < nBoundaries; boundaryI++)
    {
    const int startFace = this->BoundaryDict[boundaryI].startFace;
    const int endFace = startFace + this->BoundaryDict[boundaryI].nFaces;
    int nPoints = 0;
    for(int j = startFace; j < endFace; j++)
      {
      const int nFacePoints = facesPoints->size(j);
      nPoints += nFacePoints;
      if(nFacePoints > maxNFacePoints)
        {
        maxNFacePoints = nFacePoints;
        }
      }
    nBoundaryPointsList->SetValue(boundaryI, nPoints);
    }

  // aloocate array for converting int vector to vtkIdType List:
  // workaround for 64bit machines
  vtkIdList *facePointsVtkId = vtkIdList::New();
  facePointsVtkId->SetNumberOfIds(maxNFacePoints);

  // create initial internal point list: set all points to -1
  if(this->Parent->GetCreateCellToPoint())
    {
    this->InternalPoints = vtkIntArray::New();
    this->InternalPoints->SetNumberOfValues(this->NumPoints);
    for(int pointI = 0; pointI < this->NumPoints; pointI++)
      {
      this->InternalPoints->SetValue(pointI, -1);
      }

    // mark boundary points as 0
    for(int boundaryI = 0; boundaryI < nBoundaries; boundaryI++)
      {
      const vtkFoamBoundaryEntry &beI = this->BoundaryDict[boundaryI];
      if(beI.boundaryType == vtkFoamBoundaryEntry::PHYSICAL
        || beI.boundaryType == vtkFoamBoundaryEntry::PROCESSOR
	)
        {
        const int startFace = beI.startFace;
        const int endFace = startFace + beI.nFaces;

        for(int j = startFace; j < endFace; j++)
          {
          const int *facePoints = facesPoints->operator[](j);
          const int nFacePoints = facesPoints->size(j);
          for(int k = 0; k < nFacePoints; k++)
            {
            this->InternalPoints->SetValue(facePoints[k], 0);
            }
          }
        }
      }
    }

  int nAllBoundaryPoints = 0;
  vtkstd::vector<vtkstd::vector<int> > procCellList;
  vtkIntArray *pointTypes = NULL;

  if(this->Parent->GetCreateCellToPoint())
    {
    // create global to AllBounaries point map
    for(int pointI = 0; pointI < this->NumPoints; pointI++)
      {
      if(this->InternalPoints->GetValue(pointI) == 0)
        {
        this->InternalPoints->SetValue(pointI, nAllBoundaryPoints);
        nAllBoundaryPoints++;
        }
      }

    if(this->ProcessorName != "")
      {
      // initialize physical-processor boundary shared point list
      procCellList.resize(nAllBoundaryPoints);
      pointTypes = vtkIntArray::New();
      pointTypes->SetNumberOfTuples(nAllBoundaryPoints);
      for(int pointI = 0; pointI < nAllBoundaryPoints; pointI++)
        {
        pointTypes->SetValue(pointI, 0);
        }
      }
    }

  for(int boundaryI = 0; boundaryI < nBoundaries; boundaryI++)
    {
    const vtkFoamBoundaryEntry &beI = this->BoundaryDict[boundaryI];
    const int nFaces = beI.nFaces;
    const int startFace = beI.startFace;
    const int endFace = startFace + nFaces;

    if(this->Parent->GetCreateCellToPoint()
      && (beI.boundaryType == vtkFoamBoundaryEntry::PHYSICAL
      || beI.boundaryType == vtkFoamBoundaryEntry::PROCESSOR
      ))
      {
      // add faces to AllBoundaries
      this->InsertFacesToGrid(this->AllBoundaries, facesPoints, startFace,
        endFace, this->InternalPoints, facePointsVtkId, NULL, false);

      if(this->ProcessorName != "")
        {
        // mark belonging boundary types and, if PROCESSOR, cell numbers
        const int abStartFace = beI.allBoundariesStartFace;
        const int abEndFace = abStartFace + beI.nFaces;
        for(int faceI = abStartFace; faceI < abEndFace; faceI++)
          {
          vtkIdType nPoints;
          vtkIdType *points;
          this->AllBoundaries->GetCellPoints(faceI, nPoints, points);
          if(beI.boundaryType == vtkFoamBoundaryEntry::PHYSICAL)
            {
            for(int pointI = 0; pointI < nPoints; pointI++)
              {
              *pointTypes->GetPointer(points[pointI])
                |= vtkFoamBoundaryEntry::PHYSICAL;
              }
            }
          else // PROCESSOR
            {
            for(int pointI = 0; pointI < nPoints; pointI++)
              {
              const int pointJ = points[pointI];
              *pointTypes->GetPointer(pointJ)
                |= vtkFoamBoundaryEntry::PROCESSOR;
              procCellList[pointJ].push_back(faceI);
              }
            }
          }
	}
      }

    // skip below if inactive
    if(!beI.isActive)
      {
      continue;
      }

    // create the mesh
    const unsigned int activeBoundaryI = boundaryMesh->GetNumberOfBlocks();
    vtkUnstructuredGrid *bm = vtkUnstructuredGrid::New();
    boundaryMesh->SetBlock(activeBoundaryI, bm);

    // set the name of boundary
    this->SetBlockName(boundaryMesh, activeBoundaryI, beI.boundaryName.c_str());

    bm->Allocate(nFaces);
    const int nBoundaryPoints = nBoundaryPointsList->GetValue(boundaryI);

    // create global to boundary-local point map and boundary points
    vtkIntArray *boundaryPointList = vtkIntArray::New();
    boundaryPointList->SetNumberOfValues(nBoundaryPoints);
    int pointI = 0;
    for(int j = startFace; j < endFace; j++)
      {
      const int *facePoints = facesPoints->operator[](j);
      int nFacePoints = facesPoints->size(j);
      for(int k = 0; k < nFacePoints; k++)
        {
        boundaryPointList->SetValue(pointI, facePoints[k]);
        pointI++;
        }
      }
    vtkSortDataArray::Sort(boundaryPointList);
    this->BoundaryPointMap->push_back(vtkIntArray::New());
    vtkIntArray& bpMap = *this->BoundaryPointMap->back();
    vtkFloatArray *boundaryPointArray = vtkFloatArray::New();
    boundaryPointArray->SetNumberOfComponents(3);
    int oldPointJ = -1;
    for(int j = 0; j < nBoundaryPoints; j++)
      {
      const int pointJ = boundaryPointList->GetValue(j);
      if(pointJ != oldPointJ)
        {
        oldPointJ = pointJ;
        boundaryPointArray->InsertNextTuple(pointArray->GetPointer(3 * pointJ));
        bpMap.InsertNextValue(pointJ);
        }
      }
    boundaryPointArray->Squeeze();
    bpMap.Squeeze();
    boundaryPointList->Delete();
    vtkPoints *boundaryPoints = vtkPoints::New();
    boundaryPoints->SetData(boundaryPointArray);
    boundaryPointArray->Delete();

    // set points for boundary
    bm->SetPoints(boundaryPoints);
    boundaryPoints->Delete();

    // insert faces to boundary mesh
    this->InsertFacesToGrid(bm, facesPoints, startFace, endFace, &bpMap,
      facePointsVtkId, NULL, true);
    bm->Delete();
    bpMap.ClearLookup();
    }

  nBoundaryPointsList->Delete();
  facePointsVtkId->Delete();

  if(this->Parent->GetCreateCellToPoint())
    {
    this->AllBoundaries->Squeeze();
    this->AllBoundariesPointMap = vtkIntArray::New();
    vtkIntArray &abpMap = *this->AllBoundariesPointMap;
    abpMap.SetNumberOfValues(nAllBoundaryPoints);

    // create lists of internal points and AllBoundaries points
    int nInternalPoints = 0;
    for(int pointI = 0, allBoundaryPointI = 0; pointI < this->NumPoints;
      pointI++)
      {
      const int globalPointId = this->InternalPoints->GetValue(pointI);
      if(globalPointId == -1)
        {
        this->InternalPoints->SetValue(nInternalPoints, pointI);
        nInternalPoints++;
        }
      else
        {
        abpMap.SetValue(allBoundaryPointI, pointI);
        allBoundaryPointI++;
        }
      }
    // shrink to the number of internal points
    if(nInternalPoints > 0)
      {
      this->InternalPoints->Resize(nInternalPoints);
      }
    else
      {
      this->InternalPoints->Delete();
      this->InternalPoints = NULL;
      }

    // set dummy vtkPoints to tell the grid the number of points
    // (otherwise GetPointCells will crash)
    vtkPoints *allBoundaryPoints = vtkPoints::New();
    allBoundaryPoints->SetNumberOfPoints(abpMap.GetNumberOfTuples());
    this->AllBoundaries->SetPoints(allBoundaryPoints);
    allBoundaryPoints->Delete();

    if(this->ProcessorName != "")
      {
      // remove links to processor boundary faces from point-to-cell
      // links of physical-processor shared points to avoid cracky seams
      // on fixedValue-type boundaries which are noticeable when all the
      // decomposed meshes are appended
      this->AllBoundaries->BuildLinks();
      for(int pointI = 0; pointI < nAllBoundaryPoints; pointI++)
        {
        if(pointTypes->GetValue(pointI)
          == (vtkFoamBoundaryEntry::PHYSICAL | vtkFoamBoundaryEntry::PROCESSOR))
          {
          const vtkstd::vector<int> &procCells = procCellList[pointI];
          for(size_t cellI = 0; cellI < procCellList[pointI].size(); cellI++)
            {
            this->AllBoundaries->RemoveReferenceToCell(pointI,
              procCells[cellI]);
            }
          // omit reclaiming memory as the possibly recovered size should
          // not typically be so large
          }
        }
      pointTypes->Delete();
      }
    }

  return boundaryMesh;
}

//-----------------------------------------------------------------------------
// truncate face owner to have only boundary face info
void vtkOpenFOAMReaderPrivate::TruncateFaceOwner()
{
  const int boundaryStartFace = this->BoundaryDict[0].startFace;
  // all the boundary faces
  const int nBoundaryFaces
    = this->FaceOwner->GetNumberOfTuples() - boundaryStartFace;
  memmove(this->FaceOwner->GetPointer(0),
    this->FaceOwner->GetPointer(boundaryStartFace),
    sizeof(int) * nBoundaryFaces);
  this->FaceOwner->Resize(nBoundaryFaces);
}

//-----------------------------------------------------------------------------
// move polyhedral cell centroids
vtkPoints *vtkOpenFOAMReaderPrivate::MoveInternalMesh(
  vtkUnstructuredGrid *internalMesh, vtkFloatArray *pointArray)
{
  if(this->Parent->GetDecomposePolyhedra())
    {
    const int nAdditionalCells = this->AdditionalCellPoints->size();
    pointArray->Resize(this->NumPoints + nAdditionalCells);
    for(int i = 0; i < nAdditionalCells; i++)
      {
      vtkIntArray *polyCellPoints = this->AdditionalCellPoints->operator[](i);
      float centroid[3];
      centroid[0] = centroid[1] = centroid[2] = 0.0F;
      const int nCellPoints = polyCellPoints->GetDataSize();
      for(int j = 0; j < nCellPoints; j++)
        {
        float *pointK = pointArray->GetPointer(3 * polyCellPoints->GetValue(j));
        centroid[0] += pointK[0];
        centroid[1] += pointK[1];
        centroid[2] += pointK[2];
        }
      float weight
        = (nCellPoints ? 1.0F / static_cast<float>(nCellPoints) : 0.0F);
      centroid[0] *= weight;
      centroid[1] *= weight;
      centroid[2] *= weight;
      pointArray->InsertTuple(this->NumPoints + i, centroid);
      }
    }
  if(internalMesh->GetPoints()->GetNumberOfPoints()
    != pointArray->GetNumberOfTuples())
    {
    vtkErrorMacro(<< "The numbers of points for old points "
      << internalMesh->GetPoints()->GetNumberOfPoints() << " and new points"
      << pointArray->GetNumberOfTuples() << " don't match");
    return NULL;
    }

  // instantiate the points class
  vtkPoints* points = vtkPoints::New();
  points->SetData(pointArray);
  internalMesh->SetPoints(points);
  return points;
}

//-----------------------------------------------------------------------------
// move boundary points
void vtkOpenFOAMReaderPrivate::MoveBoundaryMesh(
  vtkMultiBlockDataSet *boundaryMesh, vtkFloatArray *pointArray)
{
  for(size_t boundaryI = 0, activeBoundaryI = 0;
    boundaryI < this->BoundaryDict.size(); boundaryI++)
    {
    if(this->BoundaryDict[boundaryI].isActive)
      {
      vtkIntArray *bpMap = this->BoundaryPointMap->operator[](activeBoundaryI);
      const int nBoundaryPoints = bpMap->GetNumberOfTuples();
      vtkFloatArray *boundaryPointArray = vtkFloatArray::New();
      boundaryPointArray->SetNumberOfComponents(3);
      boundaryPointArray->SetNumberOfTuples(nBoundaryPoints);
      for(int pointI = 0; pointI < nBoundaryPoints; pointI++)
        {
        boundaryPointArray->SetTuple(pointI, bpMap->GetValue(pointI),
          pointArray);
        }
      vtkPoints *boundaryPoints = vtkPoints::New();
      boundaryPoints->SetData(boundaryPointArray);
      boundaryPointArray->Delete();
      vtkUnstructuredGrid::SafeDownCast(boundaryMesh->GetBlock(activeBoundaryI))
        ->SetPoints(boundaryPoints);
      boundaryPoints->Delete();
      activeBoundaryI++;
      }
    }
}

//-----------------------------------------------------------------------------
// as of now the function does not do interpolation, but do just averaging.
void vtkOpenFOAMReaderPrivate:: InterpolateCellToPoint(vtkFloatArray *pData,
  vtkFloatArray *iData, vtkUnstructuredGrid *mesh, vtkIntArray *pointList,
  const int nPoints)
{
  if(nPoints == 0)
    {
    return;
    }

  // a dummy call to let GetPointCells() build the cell links if still not built
  // (not using BuildLinks() since it always rebuild links)
  vtkIdList *pointCells = vtkIdList::New();
  mesh->GetPointCells(0, pointCells);
  pointCells->Delete();

  vtkCellLinks *cl = mesh->GetCellLinks();
  const int nComponents = iData->GetNumberOfComponents();

  if(nComponents == 1)
    {
    // a special case with the innermost componentI loop unrolled
    float *tuples = iData->GetPointer(0);
    for(int pointI = 0; pointI < nPoints; pointI++)
      {
      const int pI = (pointList ? pointList->GetValue(pointI) : pointI);
      const vtkCellLinks::Link &l = cl->GetLink(pI);
      const int nCells = static_cast<int>(l.ncells);
      const vtkIdType *cells = l.cells;
      // use double intermediate variable for precision
      double interpolatedValue = 0.0;
      for(int cellI = 0; cellI < nCells; cellI++)
        {
        interpolatedValue += tuples[cells[cellI]];
        }
      interpolatedValue
        = (nCells ? interpolatedValue / static_cast<double>(nCells) : 0.0);
      pData->SetValue(pI, interpolatedValue);
      }
    }
  else if(nComponents == 3)
    {
    // a special case with the innermost componentI loop unrolled
    float *pDataPtr = pData->GetPointer(0);
    for(int pointI = 0; pointI < nPoints; pointI++)
      {
      const int pI = (pointList ? pointList->GetValue(pointI) : pointI);
      const vtkCellLinks::Link &l = cl->GetLink(pI);
      const int nCells = static_cast<int>(l.ncells);
      const vtkIdType *cells = l.cells;
      // use double intermediate variables for precision
      const double weight = (nCells ? 1.0 / static_cast<double>(nCells) : 0.0);
      double summedValue0 = 0.0, summedValue1 = 0.0, summedValue2 = 0.0;

      // hand unrolling
      for(int cellI = 0; cellI < nCells; cellI++)
        {
        const float *tuple = iData->GetPointer(3 * cells[cellI]);
        summedValue0 += tuple[0];
        summedValue1 += tuple[1];
        summedValue2 += tuple[2];
        }

      float *interpolatedValue = &pDataPtr[3 * pI];
      interpolatedValue[0] = weight * summedValue0;
      interpolatedValue[1] = weight * summedValue1;
      interpolatedValue[2] = weight * summedValue2;
      }
    }
  else
    {
    float *pDataPtr = pData->GetPointer(0);
    for(int pointI = 0; pointI < nPoints; pointI++)
      {
      const int pI = (pointList ? pointList->GetValue(pointI) : pointI);
      const vtkCellLinks::Link &l = cl->GetLink(pI);
      const int nCells = static_cast<int>(l.ncells);
      const vtkIdType *cells = l.cells;
      // use double intermediate variables for precision
      const double weight = (nCells ? 1.0 / static_cast<double>(nCells) : 0.0);
      float *interpolatedValue = &pDataPtr[nComponents * pI];
      // a bit strange loop order but this works fastest
      for(int componentI = 0; componentI < nComponents; componentI++)
        {
        const float *tuple = iData->GetPointer(componentI);
        double summedValue = 0.0;
        for(int cellI = 0; cellI < nCells; cellI++)
          {
          summedValue += tuple[nComponents * cells[cellI]];
          }
        interpolatedValue[componentI] = weight * summedValue;
        }
      }
    }
}

//-----------------------------------------------------------------------------
bool vtkOpenFOAMReaderPrivate::ReadFieldFile(vtkFoamIOobject *ioPtr,
  vtkFoamDict *dictPtr, const vtkStdString &varName,
  vtkDataArraySelection *selection)
{
  const vtkStdString varPath(this->CurrentTimeRegionPath() + "/" + varName);

  // open the file
  vtkFoamIOobject &io = *ioPtr;
  if(!io.open(varPath))
    {
    vtkErrorMacro(<<"Error opening " << io.fileName().c_str() << ": "
      << io.error().c_str());
    return false;
    }

  // if the variable is disabled on selection panel then skip it
  if(selection->ArrayExists(io.objectName().c_str())
    && !selection->ArrayIsEnabled(io.objectName().c_str()))
    {
    return false;
    }

  // read the field file into dictionary
  vtkFoamDict &dict = *dictPtr;
  if(!dict.read(io))
    {
    vtkErrorMacro(<<"Error reading line " << io.lineNumber()
      << " of " << io.fileName().c_str() << ": " << io.error().c_str());
    return false;
    }

  if(dict.type() != vtkFoamToken::DICTIONARY)
    {
    vtkErrorMacro(<<"File " << io.fileName().c_str()
      << "is not valid as a field file");
    return false;
    }
  return true;
}

//-----------------------------------------------------------------------------
vtkFloatArray *vtkOpenFOAMReaderPrivate::FillField(vtkFoamEntry *entryPtr,
  int nElements, vtkFoamIOobject *ioPtr, const vtkStdString &fieldType)
{
  vtkFloatArray *data;
  vtkFoamEntry &entry = *entryPtr;
  const vtkStdString &className = ioPtr->className();

  // "uniformValue" keyword is for uniformFixedValue B.C.
  if(entry.firstValue().isUniform() || entry.keyword() == "uniformValue")
    {
    if(entry.firstValue().type() == vtkFoamToken::SCALAR
      || entry.firstValue().type() == vtkFoamToken::LABEL)
      {
      const float num = entry.toFloat();
      data = vtkFloatArray::New();
      data->SetNumberOfValues(nElements);
      for(int i = 0; i < nElements; i++)
        {
        data->SetValue(i, num);
        }
      }
    else
      {
      float tupleBuffer[9], *tuple;
      int nComponents;
      // have to determine the type of vector
      if(entry.firstValue().type() == vtkFoamToken::LABELLIST)
        {
        vtkIntArray &ll = entry.labelList();
        nComponents = ll.GetNumberOfTuples();
        for(int componentI = 0; componentI < nComponents; componentI++)
          {
          tupleBuffer[componentI] = static_cast<float>(ll.GetValue(componentI));
          }
        tuple = tupleBuffer;
        }
      else if(entry.firstValue().type() == vtkFoamToken::SCALARLIST)
        {
        vtkFloatArray& sl = entry.scalarList();
        nComponents = sl.GetSize();
        tuple =  sl.GetPointer(0);
        }
      else
        {
        vtkErrorMacro(<<"Wrong list type for uniform field");
        return NULL;
        }

      if((fieldType == "SphericalTensorField" && nComponents == 1)
        || (fieldType == "VectorField" && nComponents == 3)
        || (fieldType == "SymmTensorField" && nComponents == 6)
        || (fieldType == "TensorField" && nComponents == 9))
        {
        data = vtkFloatArray::New();
        data->SetNumberOfComponents(nComponents);
        data->SetNumberOfTuples(nElements);
#if vtksys_DATE_STAMP_FULL >= 20080620
	// swap the components of symmTensor to match the component
	// names in paraview
	if(nComponents == 6)
	  {
	  const float symxy = tuple[1], symxz = tuple[2], symyy = tuple[3];
	  const float symyz = tuple[4], symzz = tuple[5];
	  tuple[1] = symyy;
	  tuple[2] = symzz;
	  tuple[3] = symxy;
	  tuple[4] = symyz;
	  tuple[5] = symxz;
	  }
#endif
        for(int i = 0; i < nElements; i++)
          {
          data->SetTuple(i, tuple);
          }
        }
      else
        {
        vtkErrorMacro(<< "Number of components and field class doesn't match "
          << "for " << ioPtr->objectName() << ". class = " << className
          << ", nComponents = " << nComponents);
        return NULL;
        }
      }
    }
  else // nonuniform
    {
    if((fieldType == "ScalarField"
      && entry.firstValue().type() == vtkFoamToken::SCALARLIST)
      || ((fieldType == "VectorField" || fieldType == "SphericalTensorField"
      || fieldType == "SymmTensorField" || fieldType == "TensorField")
      && entry.firstValue().type() == vtkFoamToken::VECTORLIST))
      {
      const int nTuples = entry.scalarList().GetNumberOfTuples();
      if(nTuples != nElements)
        {
        vtkErrorMacro(<<"Number of cells/points in mesh and field don't match: "
          << "mesh = " << nElements << ", field = " << nTuples);
        return NULL;
        }
      data = static_cast<vtkFloatArray *>(entry.ptr());
#if vtksys_DATE_STAMP_FULL >= 20080620
      // swap the components of symmTensor to match the component
      // names in paraview
      const int nComponents = data->GetNumberOfComponents();
      if(nComponents == 6)
	{
	for(int tupleI = 0; tupleI < nTuples; tupleI++)
	  {
	  float *tuple = data->GetPointer(nComponents * tupleI);
	  const float symxy = tuple[1], symxz = tuple[2], symyy = tuple[3];
	  const float symyz = tuple[4], symzz = tuple[5];
	  tuple[1] = symyy;
	  tuple[2] = symzz;
	  tuple[3] = symxy;
	  tuple[4] = symyz;
	  tuple[5] = symxz;
	  }
	}
#endif
      }
    else if(entry.firstValue().type() == vtkFoamToken::EMPTYLIST &&
      nElements <= 0)
      {
      data = vtkFloatArray::New();
      // set the number of components as appropriate if the list is empty
      if(fieldType == "ScalarField" || fieldType == "SphericalTensorField")
        {
        data->SetNumberOfComponents(1);
        }
      else if(fieldType == "VectorField")
        {
        data->SetNumberOfComponents(3);
        }
      else if(fieldType == "SymmTensorField")
        {
        data->SetNumberOfComponents(6);
        }
      else if(fieldType == "TensorField")
        {
        data->SetNumberOfComponents(9);
        }
      }
    else
      {
      vtkErrorMacro(<< ioPtr->objectName().c_str() << " is not a valid "
        << ioPtr->className().c_str());
      return NULL;
      }
    }
  return data;
}

//-----------------------------------------------------------------------------
// convert OpenFOAM's dimension array representation to string
void vtkOpenFOAMReaderPrivate::ConstructDimensions(vtkStdString *dimString,
  vtkFoamDict *dictPtr)
{
  if(!this->Parent->GetAddDimensionsToArrayNames())
    {
    return;
    }
  vtkFoamEntry &dimEntry = dictPtr->lookup("dimensions");
  if(dimEntry.found()
    && dimEntry.firstValue().type() == vtkFoamToken::LABELLIST)
    {
    vtkIntArray &dims = dimEntry.labelList();
    if(dims.GetNumberOfTuples() == 7)
      {
      int dimSet[7];
      for(int dimI = 0; dimI < 7; dimI++)
        {
        dimSet[dimI] = dims.GetValue(dimI);
        }
      static const char *units[7] = { "kg", "m", "s", "K", "mol", "A", "cd" };
      vtksys_ios::ostringstream posDim, negDim;
      int posSpc = 0, negSpc = 0;
      if(dimSet[0] == 1 && dimSet[1] == -1 && dimSet[2] == -2)
        {
        posDim << "Pa";
        dimSet[0] = dimSet[1] = dimSet[2] = 0;
        posSpc = 1;
        }
      for(int dimI = 0; dimI < 7; dimI++)
        {
        const int dimDim = dimSet[dimI];
        if(dimDim > 0)
          {
          if(posSpc)
            {
            posDim << " ";
            }
          posDim << units[dimI];
          if(dimDim > 1)
            {
            posDim << dimDim;
            }
          posSpc++;
          }
        else if(dimDim < 0)
          {
          if(negSpc)
            {
            negDim << " ";
            }
          negDim << units[dimI];
          if(dimDim < -1)
            {
            negDim << -dimDim;
            }
          negSpc++;
          }
        }
      *dimString += " [" + posDim.str();
      if(negSpc > 0)
        {
        if(posSpc == 0)
          {
          *dimString += "1";
          }
        if(negSpc > 1)
          {
          *dimString += "/(" + negDim.str() + ")";
          }
        else
          {
          *dimString += "/" + negDim.str();
          }
        }
      else if(posSpc == 0)
        {
        *dimString += "-";
        }
      *dimString += "]";
      }
    }
}

//-----------------------------------------------------------------------------
void vtkOpenFOAMReaderPrivate::GetVolFieldAtTimeStep(
  vtkUnstructuredGrid *internalMesh, vtkMultiBlockDataSet *boundaryMesh,
  const vtkStdString &varName)
{
  vtkFoamIOobject io(this->CasePath);
  vtkFoamDict dict;
  if(!this->ReadFieldFile(&io, &dict, varName,
    this->Parent->CellDataArraySelection))
    {
    return;
    }

  if(io.className().substr(0, 3) != "vol")
    {
    vtkErrorMacro(<< io.objectName().c_str() << " is not a volField");
    return;
    }

  vtkFoamEntry &iEntry = dict.lookup("internalField");
  if(!iEntry.found())
    {
    vtkErrorMacro(<<"internalField not found in " << io.fileName().c_str());
    return;
    }

  if(iEntry.firstValue().type() == vtkFoamToken::EMPTYLIST)
    {
    // if there's no cell there shouldn't be any boundary faces either
    if(this->NumCells > 0)
      {
      vtkErrorMacro(<<"internalField of " << io.objectName().c_str()
        << " is empty");
      }
    return;
    }

  vtkStdString fieldType = io.className().substr(3, vtkStdString::npos);
  vtkFloatArray *iData = this->FillField(&iEntry, this->NumCells, &io,
    fieldType);
  if(iData == NULL)
    {
    return;
    }

  vtkStdString dimString;
  this->ConstructDimensions(&dimString, &dict);

  vtkFloatArray *acData = NULL, *ctpData = NULL;

  if(this->Parent->GetCreateCellToPoint())
    {
    acData = vtkFloatArray::New();
    acData->SetNumberOfComponents(iData->GetNumberOfComponents());
    acData->SetNumberOfTuples(this->AllBoundaries->GetNumberOfCells());
    }

  if(iData->GetSize() > 0)
    {
    // Add field only if internal Mesh exists (skip if not selected).
    // Note we still need to read internalField even if internal mesh is
    // not selected, since boundaries without value entries may refer to
    // the internalField.
    if(internalMesh != NULL)
      {
      int nAdditionalCells = 0;
      if(this->Parent->GetDecomposePolyhedra())
        {
        // add values for decomposed cells
        nAdditionalCells = this->AdditionalCellIds->GetNumberOfTuples();
        iData->Resize(this->NumCells + nAdditionalCells);
        for(int i = 0; i < nAdditionalCells; i++)
          {
          iData->InsertTuple(this->NumCells + i,
          static_cast<vtkIdType>(this->AdditionalCellIds->GetValue(i)), iData);
          }
        }

      // set data to internal mesh
      this->AddArrayToFieldData(internalMesh->GetCellData(), iData,
	io.objectName() + dimString);

      if(this->Parent->GetCreateCellToPoint())
        {
        // Create cell-to-point interpolated data
        ctpData = vtkFloatArray::New();
        ctpData->SetNumberOfComponents(iData->GetNumberOfComponents());
        ctpData->SetNumberOfTuples(
          internalMesh->GetPoints()->GetNumberOfPoints());
        if(this->InternalPoints != NULL)
          {
          this->InterpolateCellToPoint(ctpData, iData, internalMesh,
            this->InternalPoints, this->InternalPoints->GetNumberOfTuples());
          }

        if(this->Parent->GetDecomposePolyhedra())
          {
          // assign cell values to additional points
          for(int cellI = 0, oldCellId = -1, pointI = this->NumPoints;
            cellI < nAdditionalCells; cellI++)
            {
            const int cellId = this->AdditionalCellIds->GetValue(cellI);
            if(cellId != oldCellId)
              {
              ctpData->SetTuple(pointI, cellId, iData);
              pointI++;
              oldCellId = cellId;
              }
            }
          }
        }
      }
    }
  else
    {
    // determine as there's no cells
    iData->Delete();
    if(acData != NULL)
      {
      acData->Delete();
      }
    return;
    }

  // set boundary values
  vtkFoamEntry& bEntry = dict.lookup("boundaryField");
  if(!bEntry.found())
    {
    vtkErrorMacro(<< "boundaryField not found in object " << varName.c_str()
      << " at time = " << this->TimeNames->GetValue(this->TimeStep).c_str());
    iData->Delete();
    if(acData != NULL)
      {
      acData->Delete();
      }
    if(ctpData != NULL)
      {
      ctpData->Delete();
      }
    return;
    }

  vtkstd::vector<vtkFloatArray *> vDataVector;
  for(size_t boundaryI = 0, activeBoundaryI = 0;
    boundaryI < this->BoundaryDict.size(); boundaryI++)
    {
    const vtkFoamBoundaryEntry &beI = this->BoundaryDict[boundaryI];
    const vtkStdString &boundaryNameI = beI.boundaryName;

    vtkFoamEntry& bEntryI = bEntry.dictionary().lookup(boundaryNameI);
    if(!bEntryI.found())
      {
      vtkErrorMacro(<< "boundaryField " << boundaryNameI.c_str()
        << " not found in object " << varName.c_str() << " at time = "
        << this->TimeNames->GetValue(this->TimeStep).c_str());
      iData->Delete();
      if(acData != NULL)
        {
        acData->Delete();
        }
      if(ctpData != NULL)
        {
        ctpData->Delete();
        }
      return;
      }

    if(bEntryI.firstValue().type() != vtkFoamToken::DICTIONARY)
      {
      vtkErrorMacro(<< "Type of boundaryField " << boundaryNameI.c_str()
        << " is not a subdictionary in object " << varName.c_str()
        << " at time = " << this->TimeNames->GetValue(this->TimeStep).c_str());
      iData->Delete();
      if(acData != NULL)
        {
        acData->Delete();
        }
      if(ctpData != NULL)
        {
        ctpData->Delete();
        }
      return;
      }

    const int nFaces = beI.nFaces;

    vtkFloatArray* vData = NULL;
    bool valueFound = false;
    vtkFoamEntry& vEntry = bEntryI.dictionary().lookup("value");
    if(vEntry.found()) // the boundary has a value entry
      {
      vData = this->FillField(&vEntry, nFaces, &io, fieldType);
      if(vData == NULL)
        {
        iData->Delete();
        if(acData != NULL)
          {
          acData->Delete();
          }
        if(ctpData != NULL)
          {
          ctpData->Delete();
          }
        return;
        }
      valueFound = true;
      }
    else
      {
      // uniformFixedValue B.C.
      vtkFoamEntry& ufvEntry = bEntryI.dictionary().lookup("type");
      if(ufvEntry.found())
        {
        if(ufvEntry.toString() == "uniformFixedValue")
          {
          // the boundary is of uniformFixedValue type
          vtkFoamEntry& uvEntry = bEntryI.dictionary().lookup("uniformValue");
          if(uvEntry.found()) // and has a uniformValue entry
            {
            vData = this->FillField(&uvEntry, nFaces, &io, fieldType);
            if(vData == NULL)
              {
              iData->Delete();
              if(acData != NULL)
                {
                acData->Delete();
                }
              if(ctpData != NULL)
                {
                ctpData->Delete();
                }
              return;
              }
            valueFound = true;
            }
          }
        }
      }

    const int boundaryStartFace
      = beI.startFace - this->BoundaryDict[0].startFace;

    if(!valueFound) // doesn't have a value nor uniformValue entry
      {
      // use patch-internal values as boundary values
      vData = vtkFloatArray::New();
      vData->SetNumberOfComponents(iData->GetNumberOfComponents());
      vData->SetNumberOfTuples(nFaces);
      for(int j = 0; j < nFaces; j++)
        {
        const int cellId = this->FaceOwner->GetValue(boundaryStartFace + j);
        vData->SetTuple(j, cellId, iData);
        }
      }

    if(this->Parent->GetCreateCellToPoint())
      {
      const int startFace = beI.allBoundariesStartFace;
      // if reading a processor sub-case of a decomposed case as is,
      // use the patch values of the processor patch as is
      if(beI.boundaryType == vtkFoamBoundaryEntry::PHYSICAL
        || (this->ProcessorName == ""
        && beI.boundaryType == vtkFoamBoundaryEntry::PROCESSOR))
        {
        // set the same value to AllBoundaries
        for(int faceI = 0; faceI < nFaces; faceI++)
          {
          acData->SetTuple(faceI + startFace, faceI, vData);
          }
        }
      // implies && this->ProcessorName != ""
      else if(beI.boundaryType == vtkFoamBoundaryEntry::PROCESSOR)
        {
        // average patch internal value and patch value assuming the
        // patch value to be the patchInternalField of the neighbor
        // decomposed mesh. Using double precision to avoid degrade in
        // accuracy.
        const int nComponents = vData->GetNumberOfComponents();
        for(int faceI = 0; faceI < nFaces; faceI++)
          {
          const float *vTuple = vData->GetPointer(nComponents * faceI);
          const float *iTuple = iData->GetPointer(nComponents
            * this->FaceOwner->GetValue(boundaryStartFace + faceI));
          float *acTuple
            = acData->GetPointer(nComponents * (startFace + faceI));
          for(int componentI = 0; componentI < nComponents; componentI++)
            {
            acTuple[componentI] = (static_cast<double>(vTuple[componentI])
              + static_cast<double>(iTuple[componentI])) * 0.5;
            }
          }
        }
      }

    if(beI.isActive)
      {
      vtkUnstructuredGrid *bm = vtkUnstructuredGrid::SafeDownCast(
        boundaryMesh->GetBlock(activeBoundaryI));
      this->AddArrayToFieldData(bm->GetCellData(), vData,
        io.objectName() + dimString);

      if(this->Parent->GetCreateCellToPoint())
	{
        // construct cell-to-point interpolated boundary values. This
        // is done independently from allBoundary interpolation so
        // that the interpolated values are not affected by
        // neighboring patches especially at patch edges and for
        // baffle patches
        vtkFloatArray *pData = vtkFloatArray::New();
        pData->SetNumberOfComponents(vData->GetNumberOfComponents());
        const int nPoints = bm->GetPoints()->GetNumberOfPoints();
        pData->SetNumberOfTuples(nPoints);
        this->InterpolateCellToPoint(pData, vData, bm, NULL, nPoints);
	this->AddArrayToFieldData(bm->GetPointData(), pData,
          io.objectName() + dimString);
        pData->Delete();
	}

      activeBoundaryI++;
      }
    vData->Delete();
    }
  iData->Delete();

  if(this->Parent->GetCreateCellToPoint())
    {
    // Create cell-to-point interpolated data for all boundaries and
    // override internal values
    vtkFloatArray *bpData = vtkFloatArray::New();
    bpData->SetNumberOfComponents(acData->GetNumberOfComponents());
    const int nPoints = this->AllBoundariesPointMap->GetNumberOfTuples();
    bpData->SetNumberOfTuples(nPoints);
    this->InterpolateCellToPoint(bpData, acData, this->AllBoundaries, NULL,
      nPoints);
    acData->Delete();

    if(ctpData != NULL)
      {
      // set cell-to-pint data for internal mesh
      for(int pointI = 0; pointI < nPoints; pointI++)
        {
        ctpData->SetTuple(this->AllBoundariesPointMap->GetValue(pointI), pointI,
          bpData);
        }
      this->AddArrayToFieldData(internalMesh->GetPointData(), ctpData,
	io.objectName() + dimString);
      ctpData->Delete();
      }

    bpData->Delete();
    }
}

//-----------------------------------------------------------------------------
// read point field at a timestep
void vtkOpenFOAMReaderPrivate::GetPointFieldAtTimeStep(
  vtkUnstructuredGrid *internalMesh, vtkMultiBlockDataSet *boundaryMesh,
  const vtkStdString &varName)
{
  vtkFoamIOobject io(this->CasePath);
  vtkFoamDict dict;
  if(!this->ReadFieldFile(&io, &dict, varName,
    this->Parent->PointDataArraySelection))
    {
    return;
    }

  if(io.className().substr(0, 5) != "point")
    {
    vtkErrorMacro(<< io.objectName().c_str() << " is not a pointField");
    return;
    }

  vtkFoamEntry &iEntry = dict.lookup("internalField");
  if(!iEntry.found())
    {
    vtkErrorMacro(<<"internalField not found in " << io.fileName().c_str());
    return;
    }

  if(iEntry.firstValue().type() == vtkFoamToken::EMPTYLIST)
    {
    // if there's no cell there shouldn't be any boundary faces either
    if(this->NumPoints > 0)
      {
      vtkErrorMacro(<<"internalField of " << io.objectName().c_str()
        << " is empty");
      }
    return;
    }

  vtkStdString fieldType = io.className().substr(5, vtkStdString::npos);
  vtkFloatArray *iData = this->FillField(&iEntry, this->NumPoints, &io,
    fieldType);
  if(iData == NULL)
    {
    return;
    }

  vtkStdString dimString;
  this->ConstructDimensions(&dimString, &dict);

  // AdditionalCellPoints is NULL if creation of InternalMesh had been skipped
  if(this->AdditionalCellPoints != NULL)
    {
    // point-to-cell interpolation to additional cell centroidal points
    // for decomposed cells
    const int nAdditionalPoints = this->AdditionalCellPoints->size();
    const int nComponents = iData->GetNumberOfComponents();
    iData->Resize(this->NumPoints + nAdditionalPoints);
    for(int i = 0; i < nAdditionalPoints; i++)
      {
      vtkIntArray *acp = this->AdditionalCellPoints->operator[](i);
      int nPoints = acp->GetDataSize();
      float interpolatedValue[9];
      for(int k = 0; k < nComponents; k++)
        {
        interpolatedValue[k] = 0.0F;
        }
      for(int j = 0; j < nPoints; j++)
        {
        const float *tuple = iData->GetPointer(nComponents * acp->GetValue(j));
        for(int k = 0; k < nComponents; k++)
          {
          interpolatedValue[k] += tuple[k];
          }
        }
      const float weight = 1.0F / static_cast<float>(nPoints);
      for(int k = 0; k < nComponents; k++)
        {
        interpolatedValue[k] *= weight;
        }
      iData->InsertTuple(this->NumPoints + i, interpolatedValue);
      }
    }

  if(iData->GetSize() > 0)
    {
    // Add field only if internal Mesh exists (skip if not selected).
    // Note we still need to read internalField even if internal mesh is
    // not selected, since boundaries without value entries may refer to
    // the internalField.
    if(internalMesh != NULL)
      {
      // set data to internal mesh
      this->AddArrayToFieldData(internalMesh->GetPointData(), iData,
	io.objectName() + dimString);
      }
    }
  else
    {
    // determine as there's no points
    iData->Delete();
    return;
    }

  // use patch-internal values as boundary values
  for(size_t boundaryI = 0, activeBoundaryI = 0;
    boundaryI < this->BoundaryDict.size(); boundaryI++)
    {
    if(this->BoundaryDict[boundaryI].isActive)
      {
      vtkFloatArray *vData = vtkFloatArray::New();
      vtkIntArray& bpMap = *this->BoundaryPointMap->operator[](activeBoundaryI);
      const int nPoints = bpMap.GetNumberOfTuples();
      vData->SetNumberOfComponents(iData->GetNumberOfComponents());
      vData->SetNumberOfTuples(nPoints);
      for(int j = 0; j < nPoints; j++)
        {
        vData->SetTuple(j, bpMap.GetValue(j), iData);
        }
      this->AddArrayToFieldData(vtkUnstructuredGrid::SafeDownCast(
        boundaryMesh->GetBlock(activeBoundaryI))->GetPointData(), vData,
        io.objectName() + dimString);
      vData->Delete();
      activeBoundaryI++;
      }
    }
  iData->Delete();
}

//-----------------------------------------------------------------------------
vtkMultiBlockDataSet* vtkOpenFOAMReaderPrivate::MakeLagrangianMesh()
{
  vtkMultiBlockDataSet *lagrangianMesh = vtkMultiBlockDataSet::New();

  for(int cloudI = 0;
    cloudI < this->Parent->LagrangianPaths->GetNumberOfTuples(); cloudI++)
    {
    const vtkStdString& pathI = this->Parent->LagrangianPaths->GetValue(cloudI);

    // still can't distinguish on patch selection panel, but can
    // distinguish the "lagrangian" reserved path component and a mesh
    // region with the same name
    vtkStdString subCloudName;
    if(pathI[0] == '/')
      {
      subCloudName = pathI.substr(1, vtkStdString::npos);
      }
    else
      {
      subCloudName = pathI;
      }
    if(this->RegionName != pathI.substr(0, pathI.find('/'))
      || !this->Parent->GetPatchArrayStatus(subCloudName.c_str()))
      {
      continue;
      }

    const vtkStdString cloudPath(
      this->CurrentTimePath() + "/" + subCloudName + "/");
    const vtkStdString positionsPath(cloudPath + "positions");

    // create an empty mesh to keep node/leaf structure of the
    // multi-block consistent even if mesh doesn't exist
    vtkPolyData *meshI = vtkPolyData::New();
    const int blockI = lagrangianMesh->GetNumberOfBlocks();
    lagrangianMesh->SetBlock(blockI, meshI);
    // extract the cloud name
    this->SetBlockName(lagrangianMesh, blockI,
      pathI.substr(pathI.rfind('/') + 1).c_str());

    vtkFoamIOobject io(this->CasePath);
    if(!(io.open(positionsPath) || io.open(positionsPath + ".gz")))
      {
      meshI->Delete();
      continue;
      }

    // tell the IO object if the file is in OF 1.3 binary
    // lagrangian/positions format
    io.setIs13Positions(this->Parent->GetPositionsIsIn13Format() != 0);

    vtkFoamEntryValue dict(NULL);
    try
      {
      dict.readNonuniformList<vtkFoamToken::VECTORLIST,
        vtkFoamEntryValue::vectorListTraits<vtkFloatArray, float, 3, true> >(
        io);
      }
    catch(vtkFoamError& e)
      {
      vtkErrorMacro(<<"Error reading line " << io.lineNumber()
        << " of " << io.fileName().c_str() << ": " << e.c_str());
      meshI->Delete();
      continue;
      }
    io.close();

    vtkFloatArray *pointArray = reinterpret_cast<vtkFloatArray *>(dict.ptr());
    const int nParticles = pointArray->GetNumberOfTuples();

    // instantiate the points class
    vtkPoints *points = vtkPoints::New();
    points->SetData(pointArray);
    pointArray->Delete();

    // create lagrangian mesh
    meshI->Allocate(nParticles);
    for(vtkIdType i = 0; i < nParticles; i++)
      {
      meshI->InsertNextCell(VTK_VERTEX, 1, &i);
      }
    meshI->SetPoints(points);
    points->Delete();

    // read lagrangian fields
    for(int fieldI = 0;
        fieldI < this->LagrangianFieldFiles->GetNumberOfValues(); fieldI++)
      {
      const vtkStdString varPath(
        cloudPath + this->LagrangianFieldFiles->GetValue(fieldI));

      vtkFoamIOobject io(this->CasePath);
      if(!io.open(varPath))
        {
        // if the field file doesn't exist we simply return without
        // issuing an error as a simple way of supporting multi-region
        // lagrangians
#if 0
        vtkErrorMacro(<<"Error opening " << io.fileName().c_str() << ": "
          << io.error().c_str());
#endif
        continue;
        }

      // if the variable is disabled on selection panel then skip it
      const vtkStdString selectionName(io.objectName());
      if(this->Parent->LagrangianDataArraySelection->ArrayExists(
        selectionName.c_str()) && !this->Parent->GetLagrangianArrayStatus(
          selectionName.c_str()))
        {
        continue;
        }

      // read the field file into dictionary
      vtkFoamEntryValue dict(NULL);
      if(!dict.readField(io))
        {
        vtkErrorMacro(<<"Error reading line " << io.lineNumber()
          << " of " << io.fileName().c_str() << ": " << io.error().c_str());
        continue;
        }

      // set lagrangian values
      if(dict.type() != vtkFoamToken::SCALARLIST
        && dict.type() != vtkFoamToken::VECTORLIST)
        {
        vtkErrorMacro(<< io.fileName().c_str()
          << ": Unsupported lagrangian field type " << io.className().c_str());
        continue;
        }

      vtkFloatArray* lData = reinterpret_cast<vtkFloatArray *>(dict.ptr());

      // GetNumberOfTuples() works for both scalar and vector
      const int nParticles = lData->GetNumberOfTuples();
      if(nParticles != meshI->GetNumberOfCells())
        {
        vtkErrorMacro(<< io.fileName().c_str()
          <<": Sizes of lagrangian mesh and field don't match: mesh = "
          << meshI->GetNumberOfCells() << ", field = " << nParticles);
        lData->Delete();
        continue;
        }

      this->AddArrayToFieldData(meshI->GetCellData(), lData, selectionName);
      if(this->Parent->GetCreateCellToPoint())
        {
        this->AddArrayToFieldData(meshI->GetPointData(), lData, selectionName);
        }
      lData->Delete();
      }
    meshI->Delete();
    }
  return lagrangianMesh;
}

//-----------------------------------------------------------------------------
// returns a dictionary of block names for a specified domain
vtkOpenFOAMReaderPrivate::vtkFoamDict* vtkOpenFOAMReaderPrivate::GatherBlocks(
  const char* typeIn, bool mustRead)
{
  vtkStdString type(typeIn);
  vtkStdString blockPath
    = this->CurrentTimeRegionMeshPath(this->PolyMeshFacesDir) + type;

  vtkFoamIOobject io(this->CasePath);
  if(!(io.open(blockPath) || io.open(blockPath + ".gz")))
    {
    if(mustRead)
      {
      vtkErrorMacro(<<"Error opening " << io.fileName().c_str() << ": "
        << io.error().c_str());
      }
    return NULL;
    }

  vtkFoamDict* dictPtr = new vtkFoamDict;
  vtkFoamDict& dict = *dictPtr;
  if(!dict.read(io))
    {
    vtkErrorMacro(<<"Error reading line " << io.lineNumber()
      << " of " << io.fileName().c_str() << ": " << io.error().c_str());
    delete dictPtr;
    return NULL;
    }
  if(dict.type() != vtkFoamToken::DICTIONARY)
    {
    vtkErrorMacro(<<"The file type of " << io.fileName().c_str()
      << " is not a dictionary");
    delete dictPtr;
    return NULL;
    }
  return dictPtr;
}

//-----------------------------------------------------------------------------
// returns a requested point zone mesh
bool vtkOpenFOAMReaderPrivate::GetPointZoneMesh(
  vtkMultiBlockDataSet *pointZoneMesh, vtkPoints *points)
{
  vtkFoamDict *pointZoneDictPtr = this->GatherBlocks("pointZones", false);

  if(pointZoneDictPtr == NULL)
    {
    // not an error
    return true;
    }

  vtkFoamDict &pointZoneDict = *pointZoneDictPtr;
  size_t nPointZones = pointZoneDict.size();

  for(size_t i = 0; i < nPointZones; i++)
    {
    // look up point labels
    vtkFoamDict &dict = pointZoneDict.entry(i).dictionary();
    vtkFoamEntry &pointLabelsEntry = dict.lookup("pointLabels");
    if(!pointLabelsEntry.found())
      {
      delete pointZoneDictPtr;
      vtkErrorMacro(<<"pointLabels not found in pointZones");
      return false;
      }

    // allocate an empty mesh if the list is empty
    if(pointLabelsEntry.firstValue().type() == vtkFoamToken::EMPTYLIST)
      {
      vtkUnstructuredGrid *pzm = vtkUnstructuredGrid::New();
      pointZoneMesh->SetBlock(i, pzm);
      pzm->Delete();
      // set name
      this->SetBlockName(pointZoneMesh, i,
        pointZoneDict.entry(i).keyword().c_str());
      continue;
      }

    if(pointLabelsEntry.firstValue().type() != vtkFoamToken::LABELLIST)
      {
      delete pointZoneDictPtr;
      vtkErrorMacro(<<"pointLabels not of type labelList: type = "
        << pointLabelsEntry.firstValue().type());
      return false;
      }

    vtkIntArray &labels = pointLabelsEntry.labelList();

    int nPoints = labels.GetNumberOfTuples();
    if(nPoints > this->NumPoints)
      {
      vtkErrorMacro(<<"The length of pointLabels " << nPoints
        << " for pointZone " << pointZoneDict.entry(i).keyword().c_str()
        << " exceeds the number of points " << this->NumPoints);
      delete pointZoneDictPtr;
      return false;
      }

    // allocate new grid: we do not use resize() beforehand since it
    // could lead to undefined pointer if we return by error
    vtkUnstructuredGrid *pzm = vtkUnstructuredGrid::New();

    // set pointZone size
    pzm->Allocate(nPoints);

    // insert points
    for(int j = 0; j < nPoints; j++)
      {
      vtkIdType pointLabel = labels.GetValue(j); // must be vtkIdType
      if(pointLabel >= this->NumPoints)
        {
        vtkWarningMacro(<<"pointLabels id " << pointLabel
          << " exceeds the number of points " << this->NumPoints);
        pzm->InsertNextCell(VTK_EMPTY_CELL, 0, &pointLabel);
        continue;
        }
      pzm->InsertNextCell(VTK_VERTEX, 1, &pointLabel);
      }
    pzm->SetPoints(points);

    pointZoneMesh->SetBlock(i, pzm);
    pzm->Delete();
    // set name
    this->SetBlockName(pointZoneMesh, i,
      pointZoneDict.entry(i).keyword().c_str());
    }

  delete pointZoneDictPtr;

  return true;
}

//-----------------------------------------------------------------------------
// returns a requested face zone mesh
bool vtkOpenFOAMReaderPrivate::GetFaceZoneMesh(
  vtkMultiBlockDataSet *faceZoneMesh, const intVectorVector *facesPoints,
  vtkPoints *points)
{
  vtkFoamDict *faceZoneDictPtr = this->GatherBlocks("faceZones", false);

  if(faceZoneDictPtr == NULL)
    {
    // not an error
    return true;
    }

  vtkFoamDict &faceZoneDict = *faceZoneDictPtr;
  size_t nFaceZones = faceZoneDict.size();

  for(size_t i = 0; i < nFaceZones; i++)
    {
    // look up face labels
    vtkFoamDict &dict = faceZoneDict.entry(i).dictionary();
    vtkFoamEntry &faceLabelsEntry = dict.lookup("faceLabels");
    if(!faceLabelsEntry.found())
      {
      delete faceZoneDictPtr;
      vtkErrorMacro(<<"faceLabels not found in faceZones");
      return false;
      }

    // allocate an empty mesh if the list is empty
    if(faceLabelsEntry.firstValue().type() == vtkFoamToken::EMPTYLIST)
      {
      vtkUnstructuredGrid *fzm = vtkUnstructuredGrid::New();
      faceZoneMesh->SetBlock(i, fzm);
      fzm->Delete();
      // set name
      this->SetBlockName(faceZoneMesh, i,
        faceZoneDict.entry(i).keyword().c_str());
      continue;
      }

    if(faceLabelsEntry.firstValue().type() != vtkFoamToken::LABELLIST)
      {
      delete faceZoneDictPtr;
      vtkErrorMacro(<<"faceLabels not of type labelList");
      return false;
      }

    vtkIntArray &labels = faceLabelsEntry.labelList();

    int nFaces = labels.GetNumberOfTuples();
    if(nFaces > this->FaceOwner->GetNumberOfTuples())
      {
      vtkErrorMacro(<<"The length of faceLabels " << nFaces
        << " for faceZone " << faceZoneDict.entry(i).keyword().c_str()
        << " exceeds the number of faces "
        << this->FaceOwner->GetNumberOfTuples());
      delete faceZoneDictPtr;
      return false;
      }

    // allocate new grid: we do not use resize() beforehand since it
    // could lead to undefined pointer if we return by error
    vtkUnstructuredGrid *fzm = vtkUnstructuredGrid::New();

    // set faceZone size
    fzm->Allocate(nFaces);

    // aloocate array for converting int vector to vtkIdType vector:
    // workaround for 64bit machines
    int maxNFacePoints = 0;
    for(int j = 0; j < nFaces; j++)
      {
      const int nFacePoints = facesPoints->size(labels.GetValue(j));
      if(nFacePoints > maxNFacePoints)
        {
        maxNFacePoints = nFacePoints;
        }
      }
    vtkIdList *facePointsVtkId = vtkIdList::New();
    facePointsVtkId->SetNumberOfIds(maxNFacePoints);

    // insert faces
    this->InsertFacesToGrid(fzm, facesPoints, 0, nFaces, NULL, facePointsVtkId,
      &labels, false);

    facePointsVtkId->Delete();
    fzm->SetPoints(points);
    faceZoneMesh->SetBlock(i, fzm);
    fzm->Delete();
    // set name
    this->SetBlockName(faceZoneMesh, i,
      faceZoneDict.entry(i).keyword().c_str());
    }

  delete faceZoneDictPtr;

  return true;
}

//-----------------------------------------------------------------------------
// returns a requested cell zone mesh
bool vtkOpenFOAMReaderPrivate::GetCellZoneMesh(
  vtkMultiBlockDataSet *cellZoneMesh, const intVectorVector *cellsFaces,
  const intVectorVector *facesPoints, vtkPoints *points)
{
  vtkFoamDict *cellZoneDictPtr = this->GatherBlocks("cellZones", false);

  if(cellZoneDictPtr == NULL)
    {
    // not an error
    return true;
    }

  vtkFoamDict &cellZoneDict = *cellZoneDictPtr;
  size_t nCellZones = cellZoneDict.size();

  for(size_t i = 0; i < nCellZones; i++)
    {
    // look up cell labels
    vtkFoamDict &dict = cellZoneDict.entry(i).dictionary();
    vtkFoamEntry &cellLabelsEntry = dict.lookup("cellLabels");
    if(!cellLabelsEntry.found())
      {
      delete cellZoneDictPtr;
      vtkErrorMacro(<<"cellLabels not found in cellZones");
      return false;
      }

    // allocate an empty mesh if the list is empty
    if(cellLabelsEntry.firstValue().type() == vtkFoamToken::EMPTYLIST)
      {
      vtkUnstructuredGrid *czm = vtkUnstructuredGrid::New();
      cellZoneMesh->SetBlock(i, czm);
      // set name
      this->SetBlockName(cellZoneMesh, i,
        cellZoneDict.entry(i).keyword().c_str());
      continue;
      }

    if(cellLabelsEntry.firstValue().type() != vtkFoamToken::LABELLIST)
      {
      delete cellZoneDictPtr;
      vtkErrorMacro(<<"cellLabels not of type labelList");
      return false;
      }

    vtkIntArray &labels = cellLabelsEntry.labelList();

    int nCells = labels.GetNumberOfTuples();
    if(nCells > this->NumCells)
      {
      vtkErrorMacro(<<"The length of cellLabels " << nCells
        << " for cellZone " << cellZoneDict.entry(i).keyword().c_str()
        << " exceeds the number of cells " << this->NumCells);
      delete cellZoneDictPtr;
      return false;
      }

    // allocate new grid: we do not use resize() beforehand since it
    // could lead to undefined pointers if we return by error
    vtkUnstructuredGrid *czm = vtkUnstructuredGrid::New();

    // set cellZone size
    czm->Allocate(nCells);

    // insert cells
    this->InsertCellsToGrid(czm, cellsFaces, facesPoints, NULL, NULL, &labels);

    // set cell zone points
    czm->SetPoints(points);

    cellZoneMesh->SetBlock(i, czm);
    czm->Delete();

    // set name
    this->SetBlockName(cellZoneMesh, i,
      cellZoneDict.entry(i).keyword().c_str());
    }

  delete cellZoneDictPtr;
  return true;
}

//-----------------------------------------------------------------------------
void vtkOpenFOAMReaderPrivate::AddArrayToFieldData(
  vtkDataSetAttributes *fieldData, vtkDataArray *array,
  const vtkStdString &arrayName)
{
  // exclude dimensional unit string if any
  const vtkStdString arrayNameString(arrayName.substr(0, arrayName.find(' ')));
  array->SetName(arrayName.c_str());

  if(array->GetNumberOfComponents() == 1 && arrayNameString == "p")
    {
    fieldData->SetScalars(array);
    }
  else if(array->GetNumberOfComponents() == 3 && arrayNameString == "U")
    {
    fieldData->SetVectors(array);
    }
  else
    {
    fieldData->AddArray(array);
    }
}

//-----------------------------------------------------------------------------
// return 0 if there's any error, 1 if success
int vtkOpenFOAMReaderPrivate::RequestData(vtkMultiBlockDataSet *output,
  bool recreateInternalMesh, bool recreateBoundaryMesh, bool updateVariables)
{
  recreateInternalMesh |= this->TimeStepOld == -1
    || this->InternalMeshSelectionStatus != this->InternalMeshSelectionStatusOld
    || this->PolyMeshFacesDir->GetValue(this->TimeStep)
    != this->PolyMeshFacesDir->GetValue(this->TimeStepOld)
    || this->FaceOwner == NULL;
  recreateBoundaryMesh |= recreateInternalMesh;
  updateVariables |= recreateBoundaryMesh
    || this->TimeStep != this->TimeStepOld;
  const bool pointsMoved = this->TimeStepOld == -1
    || this->PolyMeshPointsDir->GetValue(this->TimeStep)
    != this->PolyMeshPointsDir->GetValue(this->TimeStepOld);
  const bool moveInternalPoints = !recreateInternalMesh && pointsMoved;
  const bool moveBoundaryPoints = !recreateBoundaryMesh && pointsMoved;

  // RegionName check is added since subregions have region name prefixes
  const bool createEulerians = this->Parent->PatchDataArraySelection
    ->ArrayExists("internalMesh") || this->RegionName != "";

  // determine if we need to reconstruct meshes
  if(recreateInternalMesh)
    {
    this->ClearInternalMeshes();
    }
  if(recreateBoundaryMesh)
    {
    this->ClearBoundaryMeshes();
    }

  intVectorVector *facePoints = NULL;
  vtkStdString meshDir;
  if(createEulerians && (recreateInternalMesh || recreateBoundaryMesh))
    {
    // create paths to polyMesh files
    meshDir = this->CurrentTimeRegionMeshPath(this->PolyMeshFacesDir);

    // create the faces vector
    facePoints = this->ReadFacesFile(meshDir);
    if(facePoints == NULL)
      {
      return 0;
      }
    this->Parent->UpdateProgress(0.2);
    }

  intVectorVector *cellFaces = NULL;
  if(createEulerians && recreateInternalMesh)
    {
    // read owner/neighbor and create the FaceOwner and cellFaces vectors
    cellFaces = this->ReadOwnerNeighborFiles(meshDir, facePoints);
    if(cellFaces == NULL)
      {
      delete facePoints;
      return 0;
      }
    this->Parent->UpdateProgress(0.3);
    }

  vtkFloatArray *pointArray = NULL;
  if(createEulerians && (recreateInternalMesh || (recreateBoundaryMesh
    && !recreateInternalMesh && this->InternalMesh == NULL)
    || moveInternalPoints || moveBoundaryPoints))
    {
    // get the points
    pointArray = this->ReadPointsFile();
    if((pointArray == NULL && recreateInternalMesh)
      || (facePoints != NULL && !this->CheckFacePoints(facePoints)))
      {
      delete cellFaces;
      delete facePoints;
      return 0;
      }
    this->Parent->UpdateProgress(0.4);
    }

  // make internal mesh
  // Create Internal Mesh only if required for display
  if(createEulerians && recreateInternalMesh)
    {
    if(this->Parent->GetPatchArrayStatus(
      (this->RegionPrefix() + "internalMesh").c_str()))
      {
      this->InternalMesh = this->MakeInternalMesh(cellFaces, facePoints,
        pointArray);
      }
    // read and construct zones
    if(this->Parent->GetReadZones())
      {
      vtkPoints *points;
      if(this->InternalMesh != NULL)
        {
        points = this->InternalMesh->GetPoints();
        }
      else
        {
        points = vtkPoints::New();
        points->SetData(pointArray);
        }

      this->PointZoneMesh = vtkMultiBlockDataSet::New();
      if(!this->GetPointZoneMesh(this->PointZoneMesh, points))
        {
        this->PointZoneMesh->Delete();
        this->PointZoneMesh = NULL;
        delete cellFaces;
        delete facePoints;
        if(this->InternalMesh == NULL)
          {
          points->Delete();
          }
        pointArray->Delete();
        return 0;
        }
      if(this->PointZoneMesh->GetNumberOfBlocks() == 0)
        {
        this->PointZoneMesh->Delete();
        this->PointZoneMesh = NULL;
        }

      this->FaceZoneMesh = vtkMultiBlockDataSet::New();
      if(!this->GetFaceZoneMesh(this->FaceZoneMesh, facePoints, points))
        {
        this->FaceZoneMesh->Delete();
        this->FaceZoneMesh = NULL;
        if(this->PointZoneMesh != NULL)
          {
          this->PointZoneMesh->Delete();
          this->PointZoneMesh = NULL;
          }
        delete cellFaces;
        delete facePoints;
        if(this->InternalMesh == NULL)
          {
          points->Delete();
          }
        pointArray->Delete();
        return 0;
        }
      if(this->FaceZoneMesh->GetNumberOfBlocks() == 0)
        {
        this->FaceZoneMesh->Delete();
        this->FaceZoneMesh = NULL;
        }

      this->CellZoneMesh = vtkMultiBlockDataSet::New();
      if(!this->GetCellZoneMesh(this->CellZoneMesh, cellFaces, facePoints,
        points))
        {
        this->CellZoneMesh->Delete();
        this->CellZoneMesh = NULL;
        if(this->FaceZoneMesh != NULL)
          {
          this->FaceZoneMesh->Delete();
          this->FaceZoneMesh = NULL;
          }
        if(this->PointZoneMesh != NULL)
          {
          this->PointZoneMesh->Delete();
          this->PointZoneMesh = NULL;
          }
        delete cellFaces;
        delete facePoints;
        if(this->InternalMesh == NULL)
          {
          points->Delete();
          }
        pointArray->Delete();
        return 0;
        }
      if(this->CellZoneMesh->GetNumberOfBlocks() == 0)
        {
        this->CellZoneMesh->Delete();
        this->CellZoneMesh = NULL;
        }
      if(this->InternalMesh == NULL)
        {
        points->Delete();
        }
      }
    delete cellFaces;
    this->TruncateFaceOwner();
    }

  if(createEulerians && recreateBoundaryMesh)
    {
    vtkFloatArray *boundaryPointArray;
    if(pointArray != NULL)
      {
      boundaryPointArray = pointArray;
      }
    else
      {
      boundaryPointArray = static_cast<vtkFloatArray *>(
        this->InternalMesh->GetPoints()->GetData());
      }
    // create boundary mesh
    this->BoundaryMesh = this->MakeBoundaryMesh(facePoints, boundaryPointArray);
    if(this->BoundaryMesh == NULL)
      {
      delete facePoints;
      if(pointArray != NULL)
        {
        pointArray->Delete();
        }
      return 0;
      }
    }

  delete facePoints;

  // if only point coordinates change refresh point vector
  if(createEulerians && moveInternalPoints)
    {
    // refresh the points in each mesh
    vtkPoints *points;
    // Check if Internal Mesh exists first....
    if(this->InternalMesh != NULL)
      {
      points = this->MoveInternalMesh(this->InternalMesh, pointArray);
      if(points == NULL)
        {
        pointArray->Delete();
        return 0;
        }
      }
    else
      {
      points = vtkPoints::New();
      points->SetData(pointArray);
      }

    if(this->PointZoneMesh != NULL)
      {
      for(size_t i = 0; i < this->PointZoneMesh->GetNumberOfBlocks(); i++)
        {
        vtkUnstructuredGrid::SafeDownCast(this->PointZoneMesh->GetBlock(i))
          ->SetPoints(points);
        }
      }
    if(this->FaceZoneMesh != NULL)
      {
      for(size_t i = 0; i < this->FaceZoneMesh->GetNumberOfBlocks(); i++)
        {
        vtkUnstructuredGrid::SafeDownCast(this->FaceZoneMesh->GetBlock(i))
          ->SetPoints(points);
        }
      }
    if(this->CellZoneMesh != NULL)
      {
      for(size_t i = 0; i < this->CellZoneMesh->GetNumberOfBlocks(); i++)
        {
        vtkUnstructuredGrid::SafeDownCast(this->CellZoneMesh->GetBlock(i))
          ->SetPoints(points);
        }
      }
    points->Delete();
    }

  if(createEulerians && moveBoundaryPoints)
    {
    // Check if Boundary Mesh exists first....
    if(this->BoundaryMesh != NULL)
      {
      this->MoveBoundaryMesh(this->BoundaryMesh, pointArray);
      }
    }

  if(pointArray != NULL)
    {
    pointArray->Delete();
    }
  this->Parent->UpdateProgress(0.5);

  vtkMultiBlockDataSet *lagrangianMesh = NULL;
  if(updateVariables)
    {
    if(createEulerians)
      {
      if(!recreateInternalMesh && this->InternalMesh != NULL)
        {
        // clean up arrays of the previous timestep
        // Check if Internal Mesh Exists first...
        this->InternalMesh->GetCellData()->Initialize();
        this->InternalMesh->GetPointData()->Initialize();
        }
      // Check if Boundary Mesh Exists first...
      if(!recreateBoundaryMesh && this->BoundaryMesh != NULL)
        {
        for(size_t i = 0; i < this->BoundaryMesh->GetNumberOfBlocks(); i++)
          {
          vtkUnstructuredGrid *bm = vtkUnstructuredGrid::SafeDownCast(
            this->BoundaryMesh->GetBlock(i));
          bm->GetCellData()->Initialize();
          bm->GetPointData()->Initialize();
          }
        }
      // read field data variables into Internal/Boundary meshes
      for(int i = 0; i < (int)this->VolFieldFiles->GetNumberOfValues(); i++)
        {
        this->GetVolFieldAtTimeStep(this->InternalMesh, this->BoundaryMesh,
        this->VolFieldFiles->GetValue(i));
        this->Parent->UpdateProgress(0.5 + 0.25 * ((float)(i + 1)
        / ((float)this->VolFieldFiles->GetNumberOfValues() + 0.0001)));
        }
      for(int i = 0; i < (int)this->PointFieldFiles->GetNumberOfValues(); i++)
        {
        this->GetPointFieldAtTimeStep(this->InternalMesh, this->BoundaryMesh,
        this->PointFieldFiles->GetValue(i));
        this->Parent->UpdateProgress(0.75 + 0.125 * ((float)(i + 1)
        / ((float)this->PointFieldFiles->GetNumberOfValues() + 0.0001)));
        }
      }
    // read lagrangian mesh and fields
    lagrangianMesh = this->MakeLagrangianMesh();
    }

  // Add Internal Mesh to final output only if selected for display
  if(this->InternalMesh != NULL)
    {
    output->SetBlock(0, this->InternalMesh);
    this->SetBlockName(output, 0, "internalMesh");
    }

  // set boundary meshes/data as output
  if(this->BoundaryMesh != NULL && this->BoundaryMesh->GetNumberOfBlocks() > 0)
    {
    const unsigned int groupTypeI = output->GetNumberOfBlocks();
    output->SetBlock(groupTypeI, this->BoundaryMesh);
    this->SetBlockName(output, groupTypeI, "Patches");
    }

  // set lagrangian mesh as output
  if(lagrangianMesh != NULL)
    {
    if(lagrangianMesh->GetNumberOfBlocks() > 0)
      {
      const unsigned int groupTypeI = output->GetNumberOfBlocks();
      output->SetBlock(groupTypeI, lagrangianMesh);
      this->SetBlockName(output, groupTypeI, "Lagrangian Particles");
      }
    lagrangianMesh->Delete();
    }

  if(this->Parent->GetReadZones())
    {
    vtkMultiBlockDataSet *zones = NULL;
    // set Zone Meshes as output
    if(this->PointZoneMesh != NULL)
      {
      zones = vtkMultiBlockDataSet::New();
      const unsigned int zoneTypeI = zones->GetNumberOfBlocks();
      zones->SetBlock(zoneTypeI, this->PointZoneMesh);
      this->SetBlockName(zones, zoneTypeI, "pointZones");
      }

    if(this->FaceZoneMesh != NULL)
      {
      if(zones == NULL)
	{
	zones = vtkMultiBlockDataSet::New();
	}
      const unsigned int zoneTypeI = zones->GetNumberOfBlocks();
      zones->SetBlock(zoneTypeI, this->FaceZoneMesh);
      this->SetBlockName(zones, zoneTypeI, "faceZones");
      }

    if(this->CellZoneMesh != NULL)
      {
      if(zones == NULL)
	{
	zones = vtkMultiBlockDataSet::New();
	}
      const unsigned int zoneTypeI = zones->GetNumberOfBlocks();
      zones->SetBlock(zoneTypeI, this->CellZoneMesh);
      this->SetBlockName(zones, zoneTypeI, "cellZones");
      }
    if(zones != NULL)
      {
      const unsigned int groupTypeI = output->GetNumberOfBlocks();
      output->SetBlock(groupTypeI, zones);
      this->SetBlockName(output, groupTypeI, "Zones");
      }
    }

  if(this->Parent->GetCacheMesh())
    {
    this->TimeStepOld = this->TimeStep;
    }
  else
    {
    this->ClearMeshes();
    this->TimeStepOld = -1;
    }
  this->InternalMeshSelectionStatusOld = this->InternalMeshSelectionStatus;

  this->Parent->UpdateProgress(1.0);
  return 1;
}

//-----------------------------------------------------------------------------
// constructor
vtkOpenFOAMReader::vtkOpenFOAMReader()
{
  this->SetNumberOfInputPorts(0);

  this->Parent = this;
  // must be false to avoid reloading by vtkAppendCompositeDataLeaves::Update()
  this->Refresh = false;

  // INTIALIZE FILE NAME
  this->FileName = NULL;
  this->FileNameOld = new vtkStdString;

  // Child readers
  this->Readers = vtkCollection::New();

  // VTK CLASSES
  this->PatchDataArraySelection = vtkDataArraySelection::New();
  this->CellDataArraySelection = vtkDataArraySelection::New();
  this->PointDataArraySelection = vtkDataArraySelection::New();
  this->LagrangianDataArraySelection = vtkDataArraySelection::New();

  this->PatchSelectionMTimeOld = 0;
  this->CellSelectionMTimeOld = 0;
  this->PointSelectionMTimeOld = 0;
  this->LagrangianSelectionMTimeOld = 0;

  // for creating cell-to-point translated data
  this->CreateCellToPoint = 1;
  this->CreateCellToPointOld = 1;

  // for caching mesh
  this->CacheMesh = 1;

  // for decomposing polyhedra
  this->DecomposePolyhedra = 1;
  this->DecomposePolyhedraOld = 1;

  // for reading old binary lagrangian/positions format
  this->PositionsIsIn13Format = 0; // turned off by default
  this->PositionsIsIn13FormatOld = 0;

  // for reading zones
  this->ReadZones = 0; // turned off by default
  this->ReadZonesOld = 0;

  // determine if time directories are to be listed according to controlDict
  this->ListTimeStepsByControlDict = 0;
  this->ListTimeStepsByControlDictOld = 0;

  // add dimensions to array names
  this->AddDimensionsToArrayNames = 0;
  this->AddDimensionsToArrayNamesOld = 0;

  // Lagrangian paths
  this->LagrangianPaths = vtkStringArray::New();

  this->CurrentReaderIndex = 0;
  this->NumberOfReaders = 0;
}

//-----------------------------------------------------------------------------
// destructor
vtkOpenFOAMReader::~vtkOpenFOAMReader()
{
  this->LagrangianPaths->Delete();

  this->PatchDataArraySelection->Delete();
  this->CellDataArraySelection->Delete();
  this->PointDataArraySelection->Delete();
  this->LagrangianDataArraySelection->Delete();

  this->Readers->Delete();

  this->SetFileName(0);
  delete this->FileNameOld;
}

//-----------------------------------------------------------------------------
// CanReadFile
int vtkOpenFOAMReader::CanReadFile(const char *vtkNotUsed(fileName))
{
  return 1; // so far CanReadFile does nothing.
}

//-----------------------------------------------------------------------------
// PrintSelf
void vtkOpenFOAMReader::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
  os << indent << "File Name: "
     << (this->FileName ? this->FileName : "(none)") << "\n";
  return;
}

//-----------------------------------------------------------------------------
// selection list handlers

int vtkOpenFOAMReader::GetNumberOfSelectionArrays(vtkDataArraySelection *s)
{
  return s->GetNumberOfArrays();
}

int vtkOpenFOAMReader::GetSelectionArrayStatus(vtkDataArraySelection *s,
  const char *name)
{
  return s->ArrayIsEnabled(name);
}

void vtkOpenFOAMReader::SetSelectionArrayStatus(vtkDataArraySelection *s,
  const char* name, int status)
{
  unsigned long int mTime = s->GetMTime();
  if(status)
    {
    s->EnableArray(name);
    }
  else
    {
    s->DisableArray(name);
    }
  if(mTime != s->GetMTime()) // indicate that the pipeline needs to be updated
    {
    this->Modified();
    }
}

const char *vtkOpenFOAMReader::GetSelectionArrayName(vtkDataArraySelection *s,
  int index)
{
  return s->GetArrayName(index);
}

void vtkOpenFOAMReader::DisableAllSelectionArrays(vtkDataArraySelection *s)
{
  unsigned long int mTime = s->GetMTime();
  s->DisableAllArrays();
  if(mTime != s->GetMTime())
    {
    this->Modified();
    }
}

void vtkOpenFOAMReader::EnableAllSelectionArrays(vtkDataArraySelection *s)
{
  unsigned long int mTime = s->GetMTime();
  s->EnableAllArrays();
  if(mTime != s->GetMTime())
    {
    this->Modified();
    }
}

//-----------------------------------------------------------------------------
// RequestInformation
int vtkOpenFOAMReader::RequestInformation(vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *outputVector)
{
  if(!this->FileName || strlen(this->FileName) == 0)
    {
    vtkErrorMacro("FileName has to be specified!");
    return 0;
    }

  if(this->Parent == this && (*this->FileNameOld != vtkStdString(this->FileName)
    || this->ListTimeStepsByControlDict != this->ListTimeStepsByControlDictOld)
    || this->Refresh)
    {
    // clear selections
    this->CellDataArraySelection->RemoveAllArrays();
    this->PointDataArraySelection->RemoveAllArrays();
    this->LagrangianDataArraySelection->RemoveAllArrays();
    this->PatchDataArraySelection->RemoveAllArrays();

    this->NumberOfReaders = 0;

    if(!this->MakeInformationVector(outputVector, vtkStdString(""))
      || !this->MakeMetaDataAtTimeStep(true))
      {
      return 0;
      }
    this->Refresh = false;
    }
  return 1;
}

//-----------------------------------------------------------------------------
// RequestData
int vtkOpenFOAMReader::RequestData(vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *outputVector)
{
  vtkInformation* outInfo = outputVector->GetInformationObject(0);
  vtkMultiBlockDataSet *output = vtkMultiBlockDataSet::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  int nSteps = 0;
  double *requestedTimeValues = NULL;
  if(outInfo->Has(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEPS()))
    {
    requestedTimeValues
      = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEPS());
    nSteps = outInfo->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
#if 0
    double* steps
      = outInfo->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
#endif
    }

  if(nSteps > 0)
    {
    outInfo->Set(vtkDataObject::DATA_TIME_STEPS(), requestedTimeValues, 1);
    this->SetTimeValue(requestedTimeValues[0]);
    }

  if(this->Parent == this)
    {
    if(!this->MakeMetaDataAtTimeStep(false))
      {
      return 0;
      }
    this->CurrentReaderIndex = 0;
    }

  // compute flags
  // internal mesh selection change is detected within each reader
  const bool recreateInternalMesh = (!this->Parent->CacheMesh)
    || this->Parent->DecomposePolyhedra != this->Parent->DecomposePolyhedraOld
    || this->Parent->ReadZones != this->Parent->ReadZonesOld
    || this->Parent->ListTimeStepsByControlDict
    != this->Parent->ListTimeStepsByControlDictOld;
  const bool recreateBoundaryMesh =
    this->Parent->PatchDataArraySelection->GetMTime()
    != this->Parent->PatchSelectionMTimeOld
    || this->Parent->CreateCellToPoint != this->Parent->CreateCellToPointOld;
  const bool updateVariables =
    this->Parent->CellDataArraySelection->GetMTime()
    != this->Parent->CellSelectionMTimeOld
    || this->Parent->PointDataArraySelection->GetMTime()
    != this->Parent->PointSelectionMTimeOld
    || this->Parent->LagrangianDataArraySelection->GetMTime()
    != this->Parent->LagrangianSelectionMTimeOld
    || this->Parent->PositionsIsIn13Format
    != this->Parent->PositionsIsIn13FormatOld
    || this->Parent->AddDimensionsToArrayNames
    != this->Parent->AddDimensionsToArrayNamesOld;

  // create dataset
  int ret = 1;
  vtkOpenFOAMReaderPrivate *reader;
  // if the only region is not a subregion, omit being wrapped by a
  // multiblock dataset
  if(this->Readers->GetNumberOfItems() == 1
    && (reader = vtkOpenFOAMReaderPrivate::SafeDownCast(
    this->Readers->GetItemAsObject(0)))->GetRegionName() == "")
    {
    ret = reader->RequestData(output, recreateInternalMesh,
      recreateBoundaryMesh, updateVariables);
    this->Parent->CurrentReaderIndex++;
    }
  else
    {
    this->Readers->InitTraversal();
    while((reader = vtkOpenFOAMReaderPrivate::SafeDownCast(
      this->Readers->GetNextItemAsObject())) != NULL)
      {
      vtkMultiBlockDataSet *subOutput = vtkMultiBlockDataSet::New();
      if(reader->RequestData(subOutput, recreateInternalMesh,
        recreateBoundaryMesh, updateVariables))
        {
        vtkStdString regionName(reader->GetRegionName());
        if(regionName == "")
          {
          regionName = "defaultRegion";
          }
        const int blockI = output->GetNumberOfBlocks();
        output->SetBlock(blockI, subOutput);
        output->GetMetaData(blockI)->Set(vtkCompositeDataSet::NAME(),
        regionName.c_str());
        }
      else
        {
        ret = 0;
        }
      subOutput->Delete();
      this->Parent->CurrentReaderIndex++;
      }
    }

  if(this->Parent == this) // update only if this is the top-level reader
    {
    this->UpdateStatus();
    }

  return ret;
}

//-----------------------------------------------------------------------------
void vtkOpenFOAMReader::SetTimeInformation(vtkInformationVector *outputVector,
  vtkDoubleArray *timeValues)
{
  if(timeValues->GetNumberOfTuples() > 0)
    {
    outputVector->GetInformationObject(0)->Set(
      vtkStreamingDemandDrivenPipeline::TIME_STEPS(),
      timeValues->GetPointer(0), timeValues->GetNumberOfTuples());

    double timeRange[2];
    timeRange[0] = timeValues->GetValue(0);
    timeRange[1] = timeValues->GetValue(timeValues->GetNumberOfTuples() - 1);
    outputVector->GetInformationObject(0)->Set(
      vtkStreamingDemandDrivenPipeline::TIME_RANGE(), timeRange, 2);
    }
  else
    {
    double timeRange[2];
    timeRange[0] = timeRange[1] = 0.0;
    outputVector->GetInformationObject(0)->Set(
      vtkStreamingDemandDrivenPipeline::TIME_STEPS(), timeRange, 0);
    outputVector->GetInformationObject(0)->Set(
      vtkStreamingDemandDrivenPipeline::TIME_RANGE(), timeRange, 2);
    }
}

//-----------------------------------------------------------------------------
int vtkOpenFOAMReader::MakeInformationVector(vtkInformationVector *outputVector,
  const vtkStdString& procName)
{
  *this->FileNameOld = vtkStdString(this->FileName);

  // clear prior case information
  this->Readers->RemoveAllItems();

  // recreate case information
  vtkStdString casePath, controlDictPath;
  this->CreateCasePath(casePath, controlDictPath);
  casePath += procName + (procName == "" ? "" : "/");
  vtkOpenFOAMReaderPrivate *masterReader = vtkOpenFOAMReaderPrivate::New();
  if(!masterReader->MakeInformationVector(casePath, controlDictPath, procName,
    this->Parent))
    {
    masterReader->Delete();
    return 0;
    }

  if(masterReader->GetTimeValues()->GetNumberOfTuples() == 0)
    {
    vtkErrorMacro(<< this->FileName << " contains no timestep data.");
    masterReader->Delete();
    return 0;
    }

  this->Readers->AddItem(masterReader);

  if(outputVector != NULL)
    {
    this->SetTimeInformation(outputVector, masterReader->GetTimeValues());
    }

  // search subregions under constant subdirectory
  vtkStdString constantPath(casePath + "constant/");
  vtkDirectory *dir = vtkDirectory::New();
  if(!dir->Open(constantPath.c_str()))
    {
    vtkErrorMacro(<< "Can't open " << constantPath.c_str());
    return 0;
    }
  for(int fileI = 0; fileI < dir->GetNumberOfFiles(); fileI++)
    {
    vtkStdString subDir(dir->GetFile(fileI));
    if(subDir != "." && subDir != ".."
      && dir->FileIsDirectory(subDir.c_str()))
      {
      vtkStdString boundaryPath(constantPath + subDir + "/polyMesh/boundary");
      if(vtksys::SystemTools::FileExists(boundaryPath.c_str(), true)
        || vtksys::SystemTools::FileExists((boundaryPath + ".gz").c_str(),
        true))
        {
        vtkOpenFOAMReaderPrivate *subReader = vtkOpenFOAMReaderPrivate::New();
        subReader->SetupInformation(casePath, subDir, procName, masterReader);
        this->Readers->AddItem(subReader);
        subReader->Delete();
        }
      }
    }
  dir->Delete();
  masterReader->Delete();
  this->Parent->NumberOfReaders += this->Readers->GetNumberOfItems();

  return 1;
}

//-----------------------------------------------------------------------------
void vtkOpenFOAMReader::CreateCasePath(vtkStdString &casePath,
  vtkStdString &controlDictPath)
{
#if defined(_WIN32)
  const vtkStdString pathFindSeparator = "/\\", pathSeparator = "\\";
#else
  const vtkStdString pathFindSeparator = "/", pathSeparator = "/";
#endif
  controlDictPath = this->FileName;

  // determine the case directory and path to controlDict
  vtkStdString::size_type pos = controlDictPath.find_last_of(pathFindSeparator);
  if(pos == vtkStdString::npos)
    {
    // if there's no prepending path, prefix with the current directory
    controlDictPath = "." + pathSeparator + controlDictPath;
    pos = 1;
    }
  if(controlDictPath.substr(pos + 1, 11) == "controlDict")
    {
    // remove trailing "/controlDict*"
    casePath = controlDictPath.substr(0, pos - 1);
    if(casePath == ".")
      {
      casePath = ".." + pathSeparator;
      }
    else
      {
      pos = casePath.find_last_of(pathFindSeparator);
      if(pos == vtkStdString::npos)
        {
        casePath = "." + pathSeparator;
        }
      else
        {
        // remove trailing "system" (or any other directory name)
        casePath.erase(pos + 1); // preserve the last "/"
        }
      }
    }
  else
    {
    // if the file is named other than controlDict*, use the directory
    // containing the file as case directory
    casePath = controlDictPath.substr(0, pos + 1);
    controlDictPath = casePath + "system" + pathSeparator + "controlDict";
    }
}

//-----------------------------------------------------------------------------
void vtkOpenFOAMReader::AddSelectionNames(vtkDataArraySelection *selections,
  vtkStringArray *objects)
{
  objects->Squeeze();
  vtkSortDataArray::Sort(objects);
  for(int nameI = 0; nameI < objects->GetNumberOfValues(); nameI++)
    {
    selections->AddArray(objects->GetValue(nameI).c_str());
    }
  objects->Delete();
}

//-----------------------------------------------------------------------------
bool vtkOpenFOAMReader::SetTimeValue(const double timeValue)
{
  bool modified = false;
  vtkOpenFOAMReaderPrivate *reader;
  this->Readers->InitTraversal();
  while((reader = vtkOpenFOAMReaderPrivate::SafeDownCast(
    this->Readers->GetNextItemAsObject())) != NULL)
    {
    const unsigned long mTime = reader->GetMTime();
    reader->SetTimeValue(timeValue);
    if(reader->GetMTime() != mTime)
      {
      modified = true;
      }
    }
  return modified;
}

//-----------------------------------------------------------------------------
vtkDoubleArray *vtkOpenFOAMReader::GetTimeValues()
{
  vtkOpenFOAMReaderPrivate *reader = vtkOpenFOAMReaderPrivate::SafeDownCast(
    this->Readers->GetItemAsObject(0));
  return reader != NULL ? reader->GetTimeValues() : NULL;
}

//-----------------------------------------------------------------------------
int vtkOpenFOAMReader::MakeMetaDataAtTimeStep(const bool listNextTimeStep)
{
  vtkStringArray *cellSelectionNames = vtkStringArray::New();
  vtkStringArray *pointSelectionNames = vtkStringArray::New();
  vtkStringArray *lagrangianSelectionNames = vtkStringArray::New();
  int ret = 1;
  vtkOpenFOAMReaderPrivate *reader;
  this->Readers->InitTraversal();
  while((reader = vtkOpenFOAMReaderPrivate::SafeDownCast(
    this->Readers->GetNextItemAsObject())) != NULL)
    {
    ret *= reader->MakeMetaDataAtTimeStep(cellSelectionNames,
      pointSelectionNames, lagrangianSelectionNames, listNextTimeStep);
    }
  this->AddSelectionNames(this->Parent->CellDataArraySelection,
    cellSelectionNames);
  this->AddSelectionNames(this->Parent->PointDataArraySelection,
    pointSelectionNames);
  this->AddSelectionNames(this->Parent->LagrangianDataArraySelection,
    lagrangianSelectionNames);

  return ret;
}

//-----------------------------------------------------------------------------
void vtkOpenFOAMReader::UpdateStatus()
{
  // update selection MTimes
  this->PatchSelectionMTimeOld = this->PatchDataArraySelection->GetMTime();
  this->CellSelectionMTimeOld = this->CellDataArraySelection->GetMTime();
  this->PointSelectionMTimeOld = this->PointDataArraySelection->GetMTime();
  this->LagrangianSelectionMTimeOld
    = this->LagrangianDataArraySelection->GetMTime();
  this->CreateCellToPointOld = this->CreateCellToPoint;
  this->DecomposePolyhedraOld = this->DecomposePolyhedra;
  this->PositionsIsIn13FormatOld = this->PositionsIsIn13Format;
  this->ReadZonesOld = this->ReadZones;
  this->ListTimeStepsByControlDictOld = this->ListTimeStepsByControlDict;
  this->AddDimensionsToArrayNamesOld = this->AddDimensionsToArrayNames;
}

//-----------------------------------------------------------------------------
void vtkOpenFOAMReader::UpdateProgress(double amount)
{
  this->vtkAlgorithm::UpdateProgress((static_cast<double>(
    this->Parent->CurrentReaderIndex) + amount) / static_cast<double>(
    this->Parent->NumberOfReaders));
}

/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkNewOpenFOAMReader.cxx,v $

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
// Token-based FoamFile format lexer/parser,
// performance/stability/compatibility enhancements, gzipped file
// support, lagrangian field support, variable timestep support,
// builtin cell-to-point filter, pointField support, polyhedron
// decomposition support, OF 1.5 extended format support, multiregion
// support, old mesh format support, parallelization support for
// decomposed cases in conjunction with vtkNewPOpenFOAMReader, et. al. by
// Takuya Oshima of Niigata University, Japan (oshima@eng.niigata-u.ac.jp)
//
// * GUI Based selection of mesh regions and fields available in the case
// * Minor bug fixes / Strict memory allocation checks
// * Minor performance enhancements
// by Philippose Rajan (sarith@rocketmail.com)

// Hijack the CRC routine of zlib to omit CRC check for gzipped files
// (on OSes other than Windows where the mechanism doesn't work due
// to pre-bound DLL symbols) if set to 1, or not (set to 0). Affects
// performance by about 3% - 4%.
#define VTK_FOAMFILE_OMIT_CRCCHECK 0

// The input/output buffer sizes for zlib in bytes.
#define VTK_FOAMFILE_INBUFSIZE (16384)
#define VTK_FOAMFILE_OUTBUFSIZE (131072)
#define VTK_FOAMFILE_INCLUDE_STACK_SIZE (10)

// Avoid a locale problem where period is not interpreted as decimal
// point when ParaView is compiled with Qt 4.5 and run under certain
// locales (e.g. de_DE, fr_FR)
#define VTK_FOAMFILE_LOCALE_WORKAROUND 1

#if defined(_MSC_VER) && (_MSC_VER >= 1400)
#define _CRT_SECURE_NO_WARNINGS 1
#endif

#if VTK_FOAMFILE_OMIT_CRCCHECK
#define ZLIB_INTERNAL
#endif

// for possible future extension of linehead-aware directives
#define VTK_FOAMFILE_RECOGNIZE_LINEHEAD 0

#include "vtkNewOpenFOAMReader.h"

#if !defined(VTK_FOAMFILE_HAVE_REGEX)
#include "vtksys/RegularExpression.hxx"
#endif
#include <vtkstd/vector>
#include "vtksys/DateStamp.h"
#include "vtksys/SystemTools.hxx"
#include <vtksys/ios/sstream>
#include "vtk_zlib.h"

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
// for regcomp() / regexec() / regfree() / regerror()
#include <regex.h>
// for getuid()
#include <unistd.h>
#endif
// for fabs()
#include <math.h>
// for isalnum() / isspace() / isdigit()
#include <ctype.h>

#if VTK_FOAMFILE_OMIT_CRCCHECK
uLong ZEXPORT crc32(uLong, const Bytef *, uInt)
{ return 0; }
#endif

vtkCxxRevisionMacro(vtkNewOpenFOAMReader, "$Revision: 1.16 $");
vtkStandardNewMacro(vtkNewOpenFOAMReader);

// forward declarations
template <typename T> struct vtkNewFoamArrayVector
  : public vtkstd::vector<T *>
{
private:
  typedef vtkstd::vector<T *> Superclass;

public:
  ~vtkNewFoamArrayVector()
  {
    for(size_t arrayI = 0; arrayI < Superclass::size(); arrayI++)
      {
      if(Superclass::operator[](arrayI))
        {
        Superclass::operator[](arrayI)->Delete();
        }
      }
  }
};

typedef vtkNewFoamArrayVector<vtkIntArray> vtkNewFoamIntArrayVector;
typedef vtkNewFoamArrayVector<vtkFloatArray> vtkNewFoamFloatArrayVector;
struct vtkNewFoamIntVectorVector;

struct vtkNewFoamError;
struct vtkNewFoamToken;
struct vtkNewFoamFileStack;
struct vtkNewFoamFile;
struct vtkNewFoamIOobject;
template <typename T> struct vtkNewFoamReadValue;
struct vtkNewFoamEntryValue;
struct vtkNewFoamEntry;
struct vtkNewFoamDict;

//-----------------------------------------------------------------------------
// class vtkNewOpenFOAMReaderPrivate
// the reader core of vtkNewOpenFOAMReader
class vtkNewOpenFOAMReaderPrivate : public vtkObject
{
public:
  static vtkNewOpenFOAMReaderPrivate *New();
  vtkTypeRevisionMacro(vtkNewOpenFOAMReaderPrivate, vtkObject);
  void PrintSelf(ostream &, vtkIndent);

  vtkDoubleArray *GetTimeValues()
    {return this->TimeValues;}
  vtkStringArray *GetTimeNames()
    {return this->TimeNames;}
  vtkGetMacro(TimeStep, int);
  vtkSetMacro(TimeStep, int);
  const vtkStdString &GetRegionName() const
    {return this->RegionName;}

  // gather timestep information
  bool MakeInformationVector(const vtkStdString &, const vtkStdString &,
      const vtkStdString &, vtkNewOpenFOAMReader *);
  // read mesh/fields and create dataset
  int RequestData(vtkMultiBlockDataSet *, bool, bool, bool);
  void SetTimeValue(const double);
  int MakeMetaDataAtTimeStep(vtkStringArray *, vtkStringArray *,
      vtkStringArray *, const bool);
  void SetupInformation(const vtkStdString &, const vtkStdString &,
      const vtkStdString &, vtkNewOpenFOAMReaderPrivate *);

private:
  struct vtkNewFoamBoundaryEntry
    {
    enum bt
      {
      // each value must be in 2^n for bitwise operations
      PHYSICAL = 1, // patch, wall
      PROCESSOR = 2, // processor
      GEOMETRICAL = 0 // symmetryPlane, wedge, cyclic, empty, etc.
      };
    vtkStdString BoundaryName;
    int NFaces, StartFace, AllBoundariesStartFace;
    bool IsActive;
    bt BoundaryType;
    };

  struct vtkNewFoamBoundaryDict : public vtkstd::vector<vtkNewFoamBoundaryEntry>
    {
    // we need to keep the path to time directory where the current mesh
    // is read from, since boundaryDict may be accessed multiple times
    // at a timestep for patch selections
    vtkStdString TimeDir;
    };

  vtkNewOpenFOAMReader *Parent;

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
  vtkPolyData *AllBoundaries;
  vtkIntArray *AllBoundariesPointMap;
  vtkIntArray *InternalPoints;

  // for caching mesh
  vtkUnstructuredGrid *InternalMesh;
  vtkMultiBlockDataSet *BoundaryMesh;
  vtkNewFoamIntArrayVector *BoundaryPointMap;
  vtkNewFoamBoundaryDict BoundaryDict;
  vtkMultiBlockDataSet *PointZoneMesh;
  vtkMultiBlockDataSet *FaceZoneMesh;
  vtkMultiBlockDataSet *CellZoneMesh;
#if 0
  vtkNewFoamFloatArrayVector *ReciprocalDelta;
#endif

  // for polyhedra handling
  int NumTotalAdditionalCells;
  vtkIntArray *AdditionalCellIds;
  vtkIntArray *NumAdditionalCells;
  vtkNewFoamIntArrayVector *AdditionalCellPoints;

  // constructor and destructor are kept private
  vtkNewOpenFOAMReaderPrivate();
  ~vtkNewOpenFOAMReaderPrivate();

  // not implemented.
  vtkNewOpenFOAMReaderPrivate(const vtkNewOpenFOAMReaderPrivate &);
  void operator=(const vtkNewOpenFOAMReaderPrivate &);

  // clear mesh construction
  void ClearInternalMeshes();
  void ClearBoundaryMeshes();
  void ClearMeshes();

  vtkStdString RegionPath() const
    {return (this->RegionName == "" ? "" : "/") + this->RegionName;}
  vtkStdString TimePath(const int timeI) const
    {return this->CasePath + this->TimeNames->GetValue(timeI);}
  vtkStdString TimeRegionPath(const int timeI) const
    {return this->TimePath(timeI) + this->RegionPath();}
  vtkStdString TimeRegionMeshPath(vtkStringArray *dir, const int timeI) const
    {return this->CasePath + dir->GetValue(timeI) + this->RegionPath()
    + "/polyMesh/";}
  vtkStdString CurrentTimePath() const
    {return this->TimePath(this->TimeStep);}
  vtkStdString CurrentTimeRegionPath() const
    {return this->TimeRegionPath(this->TimeStep);}
  vtkStdString CurrentTimeRegionMeshPath(vtkStringArray *dir) const
    {return this->TimeRegionMeshPath(dir, this->TimeStep);}
  vtkStdString RegionPrefix() const
    {return this->RegionName + (this->RegionName == "" ? "" : "/");}

  // search time directories for mesh
  void AppendMeshDirToArray(vtkStringArray *, const vtkStdString &,
      const vtkStdString &, const int);
  void PopulatePolyMeshDirArrays();

  // search a time directory for field objects
  void GetFieldNames(const vtkStdString &, const bool, vtkStringArray *,
      vtkStringArray *);
  void SortFieldFiles(vtkStringArray *, vtkStringArray *, vtkStringArray *);
  void LocateLagrangianClouds(vtkStringArray *, const vtkStdString &);

  // read controlDict
  bool ListTimeDirectoriesByControlDict(vtkNewFoamDict *dict);
  bool ListTimeDirectoriesByInstances();

  // read mesh files
  vtkFloatArray* ReadPointsFile();
  vtkNewFoamIntVectorVector* ReadFacesFile (const vtkStdString &);
  vtkNewFoamIntVectorVector* ReadOwnerNeighborFiles(const vtkStdString &,
      vtkNewFoamIntVectorVector *);
  bool CheckFacePoints(vtkNewFoamIntVectorVector *);

  // create mesh
  void InsertCellsToGrid(vtkUnstructuredGrid *, const vtkNewFoamIntVectorVector *,
      const vtkNewFoamIntVectorVector *, vtkFloatArray *, vtkIdTypeArray *,
      vtkIntArray *);
  vtkUnstructuredGrid *MakeInternalMesh(const vtkNewFoamIntVectorVector *,
      const vtkNewFoamIntVectorVector *, vtkFloatArray *);
  bool InsertFacesToGrid(vtkPolyData *, const vtkNewFoamIntVectorVector *, int,
      int, vtkIntArray *, vtkIdList *, vtkIntArray *, const bool);
  template <typename T1, typename T2> bool ExtendArray(T1 *, const int);
  vtkMultiBlockDataSet* MakeBoundaryMesh(const vtkNewFoamIntVectorVector *,
      vtkFloatArray *);
  void SetBlockName(vtkMultiBlockDataSet *, unsigned int, const char *);
  void TruncateFaceOwner();
#if 0
  void CalculateReciprocalDelta(const vtkNewFoamIntVectorVector *);
#endif

  // move additional points for decomposed cells
  vtkPoints *MoveInternalMesh(vtkUnstructuredGrid *, vtkFloatArray *);
  void MoveBoundaryMesh(vtkMultiBlockDataSet *, vtkFloatArray *);

  // cell-to-point interpolator
  void InterpolateCellToPoint(vtkFloatArray *, vtkFloatArray *, vtkPointSet *,
      vtkIntArray *, const int);

  // read and create cell/point fields
  void ConstructDimensions(vtkStdString *, vtkNewFoamDict *);
  bool ReadFieldFile(vtkNewFoamIOobject *, vtkNewFoamDict *, const vtkStdString &,
      vtkDataArraySelection *);
  vtkFloatArray *FillField(vtkNewFoamEntry *, vtkNewFoamEntry *, int, vtkNewFoamIOobject *,
      const vtkStdString &, const bool);
  void GetVolFieldAtTimeStep(vtkUnstructuredGrid *, vtkMultiBlockDataSet *,
      const vtkStdString &);
  void GetPointFieldAtTimeStep(vtkUnstructuredGrid *, vtkMultiBlockDataSet *,
      const vtkStdString &);
  void AddArrayToFieldData(vtkDataSetAttributes *, vtkDataArray *,
      const vtkStdString &);

  // create lagrangian mesh/fields
  vtkMultiBlockDataSet *MakeLagrangianMesh();

  // create point/face/cell zones
  vtkNewFoamDict *GatherBlocks(const char *, const int);
  bool GetPointZoneMesh(vtkMultiBlockDataSet *, vtkPoints *);
  bool GetFaceZoneMesh(vtkMultiBlockDataSet *, const vtkNewFoamIntVectorVector *,
      vtkPoints *);
  bool GetCellZoneMesh(vtkMultiBlockDataSet *, const vtkNewFoamIntVectorVector *,
      const vtkNewFoamIntVectorVector *, vtkPoints *);
};

vtkCxxRevisionMacro(vtkNewOpenFOAMReaderPrivate, "$Revision: 1.16 $");
vtkStandardNewMacro(vtkNewOpenFOAMReaderPrivate);

//-----------------------------------------------------------------------------
// struct vtkNewFoamIntVectorVector
struct vtkNewFoamIntVectorVector
{
private:
  vtkIntArray *Indices, *Body;

  vtkNewFoamIntVectorVector();

public:
  ~vtkNewFoamIntVectorVector()
  {
    this->Indices->Delete();
    this->Body->Delete();
  }

  vtkNewFoamIntVectorVector(const vtkNewFoamIntVectorVector &ivv) :
    Indices(ivv.Indices), Body(ivv.Body)
  {
    this->Indices->Register(0); // vtkDataArrays do not have ShallowCopy
    this->Body->Register(0);
  }
  vtkNewFoamIntVectorVector(const int nElements, const int bodyLength) :
    Indices(vtkIntArray::New()), Body(vtkIntArray::New())
  {
    this->Indices->SetNumberOfValues(nElements + 1);
    this->Body->SetNumberOfValues(bodyLength);
  }

  // GetSize() returns all allocated size (Size) while GetDataSize()
  // returns actually used size (MaxId * nComponents)
  int GetBodySize() const
  {
    return this->Body->GetSize();
  }
  // note that vtkIntArray::Resize() allocates (current size + new
  // size) bytes if current size < new size
  void ResizeBody(const int bodyLength)
  {
    this->Body->Resize(bodyLength);
  }
  int *SetIndex(const int i, const int bodyI)
  {
    return this->Body->GetPointer(*this->Indices->GetPointer(i) = bodyI);
  }
  void SetValue(const int bodyI, int value)
  {
    this->Body->SetValue(bodyI, value);
  }
  const int *operator[](const int i) const
  {
    return this->Body->GetPointer(this->Indices->GetValue(i));
  }
  int GetSize(const int i) const
  {
    return this->Indices->GetValue(i + 1) - this->Indices->GetValue(i);
  }
  int GetNumberOfElements() const
  {
    return this->Indices->GetNumberOfTuples() - 1;
  }
  vtkIntArray *GetIndices()
  {
    return this->Indices;
  }
  vtkIntArray *GetBody()
  {
    return this->Body;
  }
};

//-----------------------------------------------------------------------------
// class vtkNewFoamError
// class for exception-carrying object
struct vtkNewFoamError : public vtkStdString
{
private:
  typedef vtkStdString Superclass;

public:
  vtkNewFoamError() :
    vtkStdString()
  {
  }
  vtkNewFoamError(const vtkNewFoamError& e) :
    vtkStdString(e)
  {
  }
  ~vtkNewFoamError()
  {
  }
  // a super-easy way to make use of operator<<()'s defined in
  // vtksys_ios::ostringstream class
  template <class T> vtkNewFoamError& operator<<(const T& t)
  {
    vtksys_ios::ostringstream os;
    os << t;
    this->Superclass::operator+=(os.str());
    return *this;
  }
};

//-----------------------------------------------------------------------------
// class vtkNewFoamToken
// token class which also works as container for list types
// - a word token is treated as a string token for simplicity
// - handles only atomic types. Handling of list types are left to the
//   derived classes.
struct vtkNewFoamToken
{
public:
  enum tokenType
    {
    // undefined type
    UNDEFINED,
    // atomic types
    PUNCTUATION, LABEL, SCALAR, WORD, STRING, IDENTIFIER,
    // vtkObject-derived list types
    STRINGLIST, LABELLIST, SCALARLIST, VECTORLIST,
    // original list types
    LABELLISTLIST, ENTRYVALUELIST, EMPTYLIST, DICTIONARY,
    // error state
    TOKEN_ERROR
    };

protected:
  tokenType Type;
  union
  {
    char Char;
    int Int;
    double Double;
    vtkStdString* String;
    vtkObjectBase *VtkObjectPtr;
    // vtkObject-derived list types
    vtkIntArray *LabelListPtr;
    vtkFloatArray *ScalarListPtr, *VectorListPtr;
    vtkStringArray *StringListPtr;
    // original list types
    vtkNewFoamIntVectorVector *LabelListListPtr;
    vtkstd::vector<vtkNewFoamEntryValue*> *EntryValuePtrs;
    vtkNewFoamDict *DictPtr;
  };

  void Clear()
  {
    if (this->Type == WORD || this->Type == STRING || this->Type == IDENTIFIER)
      {
      delete this->String;
      }
  }

  void AssignData(const vtkNewFoamToken& value)
  {
    switch (value.Type)
      {
      case PUNCTUATION:
        this->Char = value.Char;
        break;
      case LABEL:
        this->Int = value.Int;
        break;
      case SCALAR:
        this->Double = value.Double;
        break;
      case WORD:
      case STRING:
      case IDENTIFIER:
        this->String = new vtkStdString(*value.String);
        break;
        // required to suppress the 'enumeration value not handled' warning by
        // g++ when compiled with -Wall
      default:
        break;
      }
  }

public:
  vtkNewFoamToken() :
    Type(UNDEFINED)
  {
  }
  vtkNewFoamToken(const vtkNewFoamToken& value) :
    Type(value.Type)
  {
    this->AssignData(value);
  }
  ~vtkNewFoamToken()
  {
    this->Clear();
  }

  tokenType GetType() const
  {
    return this->Type;
  }

  template <typename T> bool Is() const;
  template <typename T> T To() const;
#if defined(_MSC_VER)
  // workaround for Win32-64ids-nmake70
  VTK_TEMPLATE_SPECIALIZE bool Is<int>() const;
  VTK_TEMPLATE_SPECIALIZE bool Is<float>() const;
  VTK_TEMPLATE_SPECIALIZE bool Is<double>() const;
  VTK_TEMPLATE_SPECIALIZE int To<int>() const;
  VTK_TEMPLATE_SPECIALIZE float To<float>() const;
  VTK_TEMPLATE_SPECIALIZE double To<double>() const;
#endif

  // workaround for SunOS-CC5.6-dbg
  int ToInt() const
  {
    return this->Int;
  }

  // workaround for SunOS-CC5.6-dbg
  float ToFloat() const
  {
    return this->Type == LABEL ? this->Int : this->Double;
  }

  bool IsWordOrString() const
  {
    return this->Type == WORD || this->Type == STRING;
  }

  const vtkStdString &ToStdString() const
  {
    return *this->String;
  }

  void SetBad()
  {
    this->Clear();
    this->Type = TOKEN_ERROR;
  }
  void SetWord(const vtkStdString &wordString)
  {
    this->operator=(wordString);
    this->Type = WORD;
  }
  void SetIdentifier(const vtkStdString& idString)
  {
    this->operator=(idString);
    this->Type = IDENTIFIER;
  }

  void operator=(const char value)
  {
    this->Clear();
    this->Type = PUNCTUATION;
    this->Char = value;
  }
  void operator=(const int value)
  {
    this->Clear();
    this->Type = LABEL;
    this->Int = value;
  }
  void operator=(const double value)
  {
    this->Clear();
    this->Type = SCALAR;
    this->Double = value;
  }
  void operator=(const char *value)
  {
    this->Clear();
    this->Type = STRING;
    this->String = new vtkStdString(value);
  }
  void operator=(const vtkStdString& value)
  {
    this->Clear();
    this->Type = STRING;
    this->String = new vtkStdString(value);
  }
  void operator=(const vtkNewFoamToken& value)
  {
    this->Clear();
    this->Type = value.Type;
    this->AssignData(value);
  }
  bool operator==(const char value) const
  {
    return this->Type == PUNCTUATION && this->Char == value;
  }
  bool operator==(const int value) const
  {
    return this->Type == LABEL && this->Int == value;
  }
  bool operator==(const vtkStdString& value) const
  {
    return this->Type == WORD && *this->String == value;
  }
  bool operator!=(const char value) const
  {
    return !this->operator==(value);
  }
  bool operator!=(const vtkStdString& value) const
  {
    return !this->operator==(value);
  }

  friend vtksys_ios::ostringstream& operator<<(vtksys_ios::ostringstream& str,
      const vtkNewFoamToken& value)
  {
    switch (value.GetType())
      {
      case TOKEN_ERROR:
        str << "badToken (an unexpected EOF?)";
        break;
      case PUNCTUATION:
        str << value.Char;
        break;
      case LABEL:
        str << value.Int;
        break;
      case SCALAR:
        str << value.Double;
        break;
      case WORD:
      case STRING:
        str << *value.String;
        break;
      case IDENTIFIER:
        str << "$" << *value.String;
        break;
      // required to suppress the 'enumeration value not handled' warning by
      // g++ when compiled with -Wall
      default:
        break;
      }
    return str;
  }
};

VTK_TEMPLATE_SPECIALIZE inline bool vtkNewFoamToken::Is<int>() const
{
  return this->Type == LABEL;
}

VTK_TEMPLATE_SPECIALIZE inline bool vtkNewFoamToken::Is<float>() const
{
  return this->Type == LABEL || this->Type == SCALAR;
}

VTK_TEMPLATE_SPECIALIZE inline bool vtkNewFoamToken::Is<double>() const
{
  return this->Type == SCALAR;
}

VTK_TEMPLATE_SPECIALIZE inline int vtkNewFoamToken::To<int>() const
{
  return this->Int;
}

VTK_TEMPLATE_SPECIALIZE inline float vtkNewFoamToken::To<float>() const
{
  return this->Type == LABEL ? this->Int : this->Double;
}

VTK_TEMPLATE_SPECIALIZE inline double vtkNewFoamToken::To<double>() const
{
  return this->Type == LABEL ? this->Int : this->Double;
}

//-----------------------------------------------------------------------------
// class vtkNewFoamFileStack
// list of variables that have to be saved when a file is included.
struct vtkNewFoamFileStack
{
protected:
  vtkStdString FileName;
  FILE *File;
  bool IsCompressed;
  z_stream Z;
  int ZStatus;
  int LineNumber;
#if VTK_FOAMFILE_RECOGNIZE_LINEHEAD
  bool WasNewline;
#endif

  // buffer pointers. using raw pointers for performance reason.
  unsigned char *Inbuf;
  unsigned char *Outbuf;
  unsigned char *BufPtr;
  unsigned char *BufEndPtr;

  vtkNewFoamFileStack() :
    FileName(), File(NULL), IsCompressed(false), ZStatus(Z_OK), LineNumber(0),
#if VTK_FOAMFILE_RECOGNIZE_LINEHEAD
        WasNewline(true),
#endif
        Inbuf(NULL), Outbuf(NULL), BufPtr(NULL), BufEndPtr(NULL)
  {
    this->Z.zalloc = Z_NULL;
    this->Z.zfree = Z_NULL;
    this->Z.opaque = Z_NULL;
  }

  void Reset()
  {
    // this->FileName = "";
    this->File = NULL;
    this->IsCompressed = false;
    // this->ZStatus = Z_OK;
    this->Z.zalloc = Z_NULL;
    this->Z.zfree = Z_NULL;
    this->Z.opaque = Z_NULL;
    // this->LineNumber = 0;
#if VTK_FOAMFILE_RECOGNIZE_LINEHEAD
    this->WasNewline = true;
#endif

    this->Inbuf = NULL;
    this->Outbuf = NULL;
    // this->BufPtr = NULL;
    // this->BufEndPtr = NULL;
  }

public:
  const vtkStdString& GetFileName() const
  {
    return this->FileName;
  }
  int GetLineNumber() const
  {
    return this->LineNumber;
  }
};

//-----------------------------------------------------------------------------
// class vtkNewFoamFile
// read and tokenize the input.
struct vtkNewFoamFile : public vtkNewFoamFileStack
{
private:
  typedef vtkNewFoamFileStack Superclass;

public:
  // #inputMode values
  enum inputModes
  {
    INPUT_MODE_PROTECT, INPUT_MODE_MERGE, INPUT_MODE_OVERWRITE, INPUT_MODE_WARN,
    INPUT_MODE_ERROR
  };

private:
  inputModes InputMode;

  // inclusion handling
  vtkNewFoamFileStack *Stack[VTK_FOAMFILE_INCLUDE_STACK_SIZE];
  int StackI;
  vtkStdString CasePath;

  // declare and define as private
  vtkNewFoamFile();
  bool InflateNext(unsigned char *buf, int requestSize);
  int NextTokenHead();
  // hacks to keep exception throwing / recursive codes out-of-line to make
  // putBack(), getc() and readExpecting() inline expandable
  void ThrowDuplicatedPutBackException();
  void ThrowUnexpectedEOFException();
  void ThrowUnexpectedNondigitCharExecption(const int c);
  void ThrowUnexpectedTokenException(const char, const int c);
  int ReadNext();

  void PutBack(const int c)
  {
    if (--this->Superclass::BufPtr < this->Superclass::Outbuf)
      {
      this->ThrowDuplicatedPutBackException();
      }
    *this->Superclass::BufPtr = c;
  }

  // get a character
  int Getc()
  {
    return this->Superclass::BufPtr == this->Superclass::BufEndPtr ? this->ReadNext()
        : *this->Superclass::BufPtr++;
  }

  vtkNewFoamError StackString()
  {
    vtksys_ios::ostringstream os;
    if (this->StackI > 0)
      {
      os << "\n included";

      for (int stackI = this->StackI - 1; stackI >= 0; stackI--)
        {
        os << " from line " << this->Stack[stackI]->GetLineNumber() << " of "
            << this->Stack[stackI]->GetFileName() << "\n";
        }
      os << ": ";
      }
    return vtkNewFoamError() << os.str();
  }

  bool CloseIncludedFile()
  {
    if (this->StackI == 0)
      {
      return false;
      }
    this->Clear();
    this->StackI--;
    // use the default bitwise assignment operator
    this->Superclass::operator=(*this->Stack[this->StackI]);
    delete this->Stack[this->StackI];
    return true;
  }

  void Clear()
  {
    if (this->Superclass::IsCompressed)
      {
      inflateEnd(&this->Superclass::Z);
      }

    delete [] this->Superclass::Inbuf;
    delete [] this->Superclass::Outbuf;
    this->Superclass::Inbuf = this->Superclass::Outbuf = NULL;

    if (this->Superclass::File)
      {
      fclose(this->Superclass::File);
      this->Superclass::File = NULL;
      }
    // don't reset the line number so that the last line number is
    // retained after close
    // lineNumber_ = 0;
  }

  const vtkStdString ExtractPath(const vtkStdString& path) const
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

  const vtkStdString ExtractCaseName() const
  {
#if defined(_WIN32)
    const vtkStdString pathFindSeparator = "/\\", pathSeparator = "\\";
#else
    const vtkStdString pathFindSeparator = "/", pathSeparator = "/";
#endif
    const vtkStdString::size_type len = this->CasePath.length();
    if (len == 0 || (len == 1
      && vtkStdString::npos != this->CasePath.find_first_of(pathFindSeparator)))
      {
      return vtkStdString();
      }
    const vtkStdString::size_type pos = this->CasePath.find_last_of(pathFindSeparator, len - 2);
    return pos == vtkStdString::npos ? this->CasePath.substr(0, len - 1)
      : this->CasePath.substr(pos + 1, len - 2 - pos);
  }

public:
  vtkNewFoamFile(const vtkStdString& casePath) :
    vtkNewFoamFileStack(), InputMode(INPUT_MODE_MERGE), StackI(0),
    CasePath(casePath)
  {
  }
  ~vtkNewFoamFile()
  {
    this->Close();
  }

  inputModes GetInputMode() const
  {
    return this->InputMode;
  }
  const vtkStdString GetCasePath() const
  {
    return this->CasePath;
  }
  const vtkStdString GetFilePath() const
  {
    return this->ExtractPath(this->FileName);
  }

  vtkStdString ExpandPath(const vtkStdString& pathIn,
      const vtkStdString& defaultPath)
  {
    vtkStdString expandedPath;
    bool isExpanded = false, wasPathSeparator = true;
    const size_t nChars = pathIn.length();
    for (size_t charI = 0; charI < nChars;)
      {
      char c = pathIn[charI];
      switch (c)
        {
        case '$': // $-variable expansion
          {
          vtkStdString variable;
          while (++charI < nChars && (isalnum(pathIn[charI]) || pathIn[charI]
              == '_'))
            {
            variable += pathIn[charI];
            }
          if (variable == "FOAM_CASE") // discard path until the variable
            {
            expandedPath = this->CasePath;
            wasPathSeparator = true;
            isExpanded = true;
            }
          else if (variable == "FOAM_CASENAME")
            {
            expandedPath += this->ExtractCaseName();
            wasPathSeparator = false;
            }
          else
            {
            const char *value = getenv(variable.c_str());
            if (value != NULL)
              {
              expandedPath += value;
              }
            const vtkStdString::size_type len = expandedPath.length();
            if (len > 0)
              {
              const char c2 = expandedPath[len - 1];
              wasPathSeparator = (c2 == '/' || c2 == '\\');
              }
            else
              {
              wasPathSeparator = false;
              }
            }
          }
          break;
        case '~': // home directory expansion
          // not using vtksys::SystemTools::ConvertToUnixSlashes() for
          // a bit better handling of "~"
          if (wasPathSeparator)
            {
            vtkStdString userName;
            while (++charI < nChars && (pathIn[charI] != '/' && pathIn[charI]
                != '\\') && pathIn[charI] != '$')
              {
              userName += pathIn[charI];
              }
            if (userName == "")
              {
              const char *homePtr = getenv("HOME");
              if (homePtr == NULL)
                {
#if defined(_WIN32) && !defined(__CYGWIN__) || defined(__LIBCATAMOUNT__)
                expandedPath = "";
#else
                const struct passwd *pwentry = getpwuid(getuid());
                if (pwentry == NULL)
                  {
                  throw this->StackString() << "Home directory path not found";
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
              expandedPath
              = this->ExtractPath(homePtr ? homePtr : "") + userName;
#else
              if (userName == "OpenFOAM")
                {
                // so far only "~/.OpenFOAM" expansion is supported
                const char *homePtr = getenv("HOME");
                if (homePtr == NULL)
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
                if (pwentry == NULL)
                  {
                  throw this->StackString() << "Home directory for user "
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
    if (isExpanded || expandedPath.substr(0, 1) == "/" || expandedPath.substr(
        0, 1) == "\\")
      {
      return expandedPath;
      }
    else
      {
      return defaultPath + expandedPath;
      }
  }

  void IncludeFile(const vtkStdString& includedFileName,
    const vtkStdString& defaultPath, const bool includeIfPresent)
  {
    if (this->StackI >= VTK_FOAMFILE_INCLUDE_STACK_SIZE)
      {
      throw this->StackString() << "Exceeded maximum #include recursions of "
      << VTK_FOAMFILE_INCLUDE_STACK_SIZE;
      }
    // use the default bitwise copy constructor
    this->Stack[this->StackI++] = new vtkNewFoamFileStack(*this);
    this->Superclass::Reset();

    if(includeIfPresent)
      {
      // the behavior of #includeIfPresent is in fact "includeIfOpenable"
      try
        {
        this->Open(this->ExpandPath(includedFileName, defaultPath));
        }
      catch(vtkNewFoamError)
        {
        this->CloseIncludedFile();
        }
      }
    else
      {
      this->Open(this->ExpandPath(includedFileName, defaultPath));
      }
  }

  // the tokenizer
  // returns true if success, false if encountered EOF
  bool Read(vtkNewFoamToken& token)
  {
    // expanded the outermost loop in nextTokenHead() for performance
    int c;
    while (isspace(c = this->Getc())) // isspace() accepts -1 as EOF
      {
      if (c == '\n')
        {
        ++this->Superclass::LineNumber;
#if VTK_FOAMFILE_RECOGNIZE_LINEHEAD
        this->Superclass::WasNewline = true;
#endif
        }
      }
    if (c == 47) // '/' == 47
      {
      this->PutBack(c);
      c = this->NextTokenHead();
      }
#if VTK_FOAMFILE_RECOGNIZE_LINEHEAD
    if(c != '#')
      {
      this->Superclass::WasNewline = false;
      }
#endif

    const int MAXLEN = 1024;
    char buf[MAXLEN + 1];
    int charI = 0;
    switch (c)
      {
      case '(':
      case ')':
        // high-priority punctuation token
        token = static_cast<char>(c);
        return true;
      case '1':
      case '2':
      case '3':
      case '4':
      case '5':
      case '6':
      case '7':
      case '8':
      case '9':
      case '0':
      case '-':
        // undetermined number token
        do
          {
          buf[charI++] = c;
          } while (isdigit(c = this->Getc()) && charI < MAXLEN);
        if (c != '.' && c != 'e' && c != 'E' && charI < MAXLEN && c != EOF)
          {
          // label token
          buf[charI] = '\0';
          token = static_cast<int>(strtol(buf, NULL, 10));
          this->PutBack(c);
          return true;
          }
        // fall through
      case '.':
        // scalar token
        if (c == '.' && charI < MAXLEN)
          {
          // read decimal fraction part
          buf[charI++] = c;
          while (isdigit(c = this->Getc()) && charI < MAXLEN)
            {
            buf[charI++] = c;
            }
          }
        if ((c == 'e' || c == 'E') && charI < MAXLEN)
          {
          // read exponent part
          buf[charI++] = c;
          if (((c = this->Getc()) == '+' || c == '-') && charI < MAXLEN)
            {
            buf[charI++] = c;
            c = this->Getc();
            }
          while (isdigit(c) && charI < MAXLEN)
            {
            buf[charI++] = c;
            c = this->Getc();
            }
          }
        if (charI == 1 && buf[0] == '-')
          {
          token = '-';
          this->PutBack(c);
          return true;
          }
        buf[charI] = '\0';
#if VTK_FOAMFILE_LOCALE_WORKAROUND
          {
          vtksys_ios::istringstream conversionStream(buf);
          double value;
          conversionStream >> value;
          token = value;
          }
#else
        token = strtod(buf, NULL);
#endif
        this->PutBack(c);
        break;
      case ';':
      case '{':
      case '}':
      case '[':
      case ']':
      case ':':
      case ',':
      case '=':
      case '+':
      case '*':
      case '/':
        // low-priority punctuation token
        token = static_cast<char>(c);
        return true;
      case '"':
        {
        // string token
        bool wasEscape = false;
        while ((c = this->Getc()) != EOF && charI < MAXLEN)
          {
          if (c == '\\' && !wasEscape)
            {
            wasEscape = true;
            continue;
            }
          else if (c == '"' && !wasEscape)
            {
            break;
            }
          else if (c == '\n')
            {
            ++this->Superclass::LineNumber;
            if (!wasEscape)
              {
              throw this->StackString()
              << "Unescaped newline in string constant";
              }
            }
          else if (wasEscape && c != '"' && c != '\n')
            {
            buf[charI++] = '\\';
            if (charI >= MAXLEN)
              {
              break;
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
        token.SetBad();
        return false;
      case '$':
        {
        vtkNewFoamToken identifierToken;
        if (!this->Read(identifierToken))
          {
          throw this->StackString() << "Unexpected EOF reading identifier";
          }
        if (identifierToken.GetType() != vtkNewFoamToken::WORD)
          {
          throw this->StackString() << "Expected a word, found "
          << identifierToken;
          }
        token.SetIdentifier(identifierToken.ToStdString());
        return true;
        }
      case '#':
        {
#if VTK_FOAMFILE_RECOGNIZE_LINEHEAD
        // placing #-directives in the middle of a line looks like
        // valid for the genuine OF 1.5 parser
        if(!this->Superclass::WasNewline)
          {
          throw this->StackString()
          << "Encountered #-directive in the middle of a line";
          }
        this->Superclass::WasNewline = false;
#endif
        // read directive
        vtkNewFoamToken directiveToken;
        if (!this->Read(directiveToken))
          {
          throw this->StackString() << "Unexpected EOF reading directive";
          }
        if (directiveToken == "include" || directiveToken == "includeIfPresent")
          {
          vtkNewFoamToken fileNameToken;
          if (!this->Read(fileNameToken))
            {
            throw this->StackString() << "Unexpected EOF reading filename";
            }
          if(fileNameToken.GetType() != vtkNewFoamToken::STRING)
            {
            throw this->StackString() << "Expected string, found " << fileNameToken;
            }
          this->IncludeFile(fileNameToken.ToStdString(), this->GetFilePath(),
            directiveToken == "includeIfPresent");
          }
        else if (directiveToken == "inputMode")
          {
          vtkNewFoamToken modeToken;
          if (!this->Read(modeToken))
            {
            throw this->StackString()
            << "Unexpected EOF reading inputMode specifier";
            }
          if (modeToken == "merge" || modeToken == "default")
            {
            this->InputMode = INPUT_MODE_MERGE;
            }
          else if (modeToken == "protect")
            {
            this->InputMode = INPUT_MODE_PROTECT;
            }
          else if (modeToken == "overwrite")
            {
            this->InputMode = INPUT_MODE_OVERWRITE;
            }
          else if (modeToken == "warn")
            {
            this->InputMode = INPUT_MODE_WARN;
            }
          else if (modeToken == "error")
            {
            this->InputMode = INPUT_MODE_ERROR;
            }
          else
            {
            throw this->StackString() << "Expected one of inputMode specifiers "
            "(merge, protect, overwrite, warn, error, default), found " << modeToken;
            }
          }
        else
          {
          throw this->StackString() << "Unsupported directive "
          << directiveToken;
          }
        return this->Read(token);
        }
      default:
        // word token
        int inBrace = 0;
        do
          {
          if (c == '(')
            {
            inBrace++;
            }
          else if (c == ')' && --inBrace == -1)
            {
            break;
            }
          buf[charI++] = c;
          // valid characters that constitutes a word
          // cf. src/OpenFOAM/primitives/strings/word/wordI.H
          } while ((c = this->Getc()) != EOF && !isspace(c) && c != '"' && c
            != '/' && c != ';' && c != '{' && c != '}' && charI < MAXLEN);
        buf[charI] = '\0';
        token.SetWord(buf);
        this->PutBack(c);
      }

    if (c == EOF)
      {
      this->ThrowUnexpectedEOFException();
      }
    if (charI == MAXLEN)
      {
      throw this->StackString() << "Exceeded maximum allowed length of "
      << MAXLEN << " chars";
      }
    return true;
  }

  void Open(const vtkStdString& fileName)
  {
    // reset line number to indicate the beginning of the file when an
    // exception is thrown
    this->Superclass::LineNumber = 0;
    this->Superclass::FileName = fileName;

    if (this->Superclass::File)
      {
      throw this->StackString() << "File already opened within this object";
      }

    if ((this->Superclass::File = fopen(this->Superclass::FileName.c_str(),
        "rb")) == NULL)
      {
      throw this->StackString() << "Can't open";
      }

    unsigned char zMagic[2];
    if (fread(zMagic, 1, 2, this->Superclass::File) == 2 && zMagic[0] == 0x1f
        && zMagic[1] == 0x8b)
      {
      // gzip-compressed format
      this->Superclass::Z.avail_in = 0;
      this->Superclass::Z.next_in = Z_NULL;
      // + 32 to automatically recognize gzip format
      if (inflateInit2(&this->Superclass::Z, 15 + 32) == Z_OK)
        {
        this->Superclass::IsCompressed = true;
        this->Superclass::Inbuf = new unsigned char[VTK_FOAMFILE_INBUFSIZE];
        }
      else
        {
        fclose(this->Superclass::File);
        this->Superclass::File = NULL;
        throw this->StackString() << "Can't init zstream "
        << (this->Superclass::Z.msg ? this->Superclass::Z.msg : "");
        }
      }
    else
      {
      // uncompressed format
      this->Superclass::IsCompressed = false;
      }
    rewind(this->Superclass::File);

    this->Superclass::ZStatus = Z_OK;
    this->Superclass::Outbuf = new unsigned char[VTK_FOAMFILE_OUTBUFSIZE + 1];
    this->Superclass::BufPtr = this->Superclass::Outbuf + 1;
    this->Superclass::BufEndPtr = this->Superclass::BufPtr;
    this->Superclass::LineNumber = 1;
  }

  void Close()
  {
    while (this->CloseIncludedFile())
      ;
    this->Clear();
  }

  // gzread with buffering handling
  int Read(unsigned char *buf, const int len)
  {
    int readlen;
    const int buflen = this->Superclass::BufEndPtr - this->Superclass::BufPtr;
    if (len > buflen)
      {
      memcpy(buf, this->Superclass::BufPtr, buflen);
      readlen = this->InflateNext(buf + buflen, len - buflen);
      if (readlen >= 0)
        {
        readlen += buflen;
        }
      else
        {
        if (buflen == 0) // return EOF
          {
          readlen = -1;
          }
        else
          {
          readlen = buflen;
          }
        }
      this->Superclass::BufPtr = this->Superclass::BufEndPtr;
      }
    else
      {
      memcpy(buf, this->Superclass::BufPtr, len);
      this->Superclass::BufPtr += len;
      readlen = len;
      }
    for (int i = 0; i < readlen; i++)
      {
      if (buf[i] == '\n')
        {
        this->Superclass::LineNumber++;
        }
      }
    return readlen;
  }

  void ReadExpecting(const char expected)
  {
    // skip prepending invalid chars
    // expanded the outermost loop in nextTokenHead() for performance
    int c;
    while (isspace(c = this->Getc())) // isspace() accepts -1 as EOF
      {
      if (c == '\n')
        {
        ++this->Superclass::LineNumber;
#if VTK_FOAMFILE_RECOGNIZE_LINEHEAD
        this->Superclass::WasNewline = true;
#endif
        }
      }
    if (c == 47) // '/' == 47
      {
      this->PutBack(c);
      c = this->NextTokenHead();
      }
    if (c != expected)
      {
      this->ThrowUnexpectedTokenException(expected, c);
      }
  }

  void ReadExpecting(const char* str)
  {
    vtkNewFoamToken t;
    if (!this->Read(t) || t != str)
      {
      throw this->StackString() << "Expected word \"" << str << "\", found "
      << t;
      }
  }

  int ReadIntValue();
  float ReadFloatValue();
};

int vtkNewFoamFile::ReadNext()
{
  if (!this->InflateNext(this->Superclass::Outbuf + 1, VTK_FOAMFILE_OUTBUFSIZE))
    {
    return this->CloseIncludedFile() ? this->Getc() : EOF;
    }
  return *this->Superclass::BufPtr++;
}

// specialized for reading an integer value.
// not using the standard strtol() for speed reason.
int vtkNewFoamFile::ReadIntValue()
{
  // skip prepending invalid chars
  // expanded the outermost loop in nextTokenHead() for performance
  int c;
  while (isspace(c = this->Getc())) // isspace() accepts -1 as EOF
    {
    if (c == '\n')
      {
      ++this->Superclass::LineNumber;
#if VTK_FOAMFILE_RECOGNIZE_LINEHEAD
      this->Superclass::WasNewline = true;
#endif
      }
    }
  if (c == 47) // '/' == 47
    {
    this->PutBack(c);
    c = this->NextTokenHead();
    }

  int nonNegative = c - 45; // '-' == 45
  if (nonNegative == 0 || c == 43) // '+' == 43
    {
    c = this->Getc();
    if (c == '\n')
      {
      ++this->Superclass::LineNumber;
#if VTK_FOAMFILE_RECOGNIZE_LINEHEAD
      this->Superclass::WasNewline = true;
#endif
      }
    }

  if (!isdigit(c)) // isdigit() accepts -1 as EOF
    {
    if (c == EOF)
      {
      this->ThrowUnexpectedEOFException();
      }
    else
      {
      this->ThrowUnexpectedNondigitCharExecption(c);
      }
    }

  int num = c - 48; // '0' == 48
  while (isdigit(c = this->Getc()))
    {
    num = 10 * num + c - 48;
    }

  if (c == EOF)
    {
    this->ThrowUnexpectedEOFException();
    }
  this->PutBack(c);

  return nonNegative ? num : -num;
}

// extreamely simplified high-performing string to floating point
// conversion code based on
// ParaView3/VTK/Utilities/vtksqlite/vtk_sqlite3.c
float vtkNewFoamFile::ReadFloatValue()
{
  // skip prepending invalid chars
  // expanded the outermost loop in nextTokenHead() for performance
  int c;
  while (isspace(c = this->Getc())) // isspace() accepts -1 as EOF
    {
    if (c == '\n')
      {
      ++this->Superclass::LineNumber;
#if VTK_FOAMFILE_RECOGNIZE_LINEHEAD
      this->Superclass::WasNewline = true;
#endif
      }
    }
  if (c == 47) // '/' == 47
    {
    this->PutBack(c);
    c = this->NextTokenHead();
    }

  // determine sign
  int nonNegative = c - 45; // '-' == 45
  if (nonNegative == 0 || c == 43) // '+' == 43
    {
    c = this->Getc();
    if (c == '\n')
      {
      ++this->Superclass::LineNumber;
#if VTK_FOAMFILE_RECOGNIZE_LINEHEAD
      this->Superclass::WasNewline = true;
#endif
      }
    }

  if (!isdigit(c) && c != 46) // '.' == 46, isdigit() accepts EOF
    {
    this->ThrowUnexpectedNondigitCharExecption(c);
    }

  // read integer part
  double num = c - 48; // '0' == 48
  while (isdigit(c = this->Getc()))
    {
    num = num * 10.0 + (c - 48);
    }

  // read decimal part
  if (c == 46) // '.'
    {
    double divisor = 1.0;

    while (isdigit(c = this->Getc()))
      {
      num = num * 10.0 + (c - 48);
      divisor *= 10.0;
      }
    num /= divisor;
    }

  // read exponent part
  if (c == 69 || c == 101) // 'E' == 69, 'e' == 101
    {
    int esign = 1;
    int eval = 0;
    double scale = 1.0;

    c = this->Getc();
    if (c == 45) // '-'
      {
      esign = -1;
      c = this->Getc();
      }
    else if (c == 43) // '+'
      {
      c = this->Getc();
      }

    while (isdigit(c))
      {
      eval = eval * 10 + (c - 48);
      c = this->Getc();
      }

    // fast exponent multiplication!
    while (eval >= 64)
      {
      scale *= 1.0e+64;
      eval -= 64;
      }
    while (eval >= 16)
      {
      scale *= 1.0e+16;
      eval -= 16;
      }
    while (eval >= 4)
      {
      scale *= 1.0e+4;
      eval -= 4;
      }
    while (eval >= 1)
      {
      scale *= 1.0e+1;
      eval -= 1;
      }

    if (esign < 0)
      {
      num /= scale;
      }
    else
      {
      num *= scale;
      }
    }

  if (c == EOF)
    {
    this->ThrowUnexpectedEOFException();
    }
  this->PutBack(c);

  return static_cast<float>(nonNegative ? num : -num);
}

// hacks to keep exception throwing code out-of-line to make
// putBack() and readExpecting() inline expandable
void vtkNewFoamFile::ThrowUnexpectedEOFException()
{
  throw this->StackString() << "Unexpected EOF";
}

void vtkNewFoamFile::ThrowUnexpectedNondigitCharExecption(const int c)
{
  throw this->StackString() << "Expected a number, found a non-digit character "
  << static_cast<char>(c);
}

void vtkNewFoamFile::ThrowUnexpectedTokenException(const char expected, const int c)
{
  vtkNewFoamError sstr;
  sstr << this->StackString() << "Expected punctuation token '" << expected
      << "', found ";
  if (c == EOF)
    {
    sstr << "EOF";
    }
  else
    {
    sstr << static_cast<char>(c);
    }
  throw sstr;
}

void vtkNewFoamFile::ThrowDuplicatedPutBackException()
{
  throw this->StackString() << "Attempted duplicated putBack()";
}

bool vtkNewFoamFile::InflateNext(unsigned char *buf,
    int requestSize)
{
  size_t size;
  if (this->Superclass::IsCompressed)
    {
    if (this->Superclass::ZStatus != Z_OK)
      {
      return false;
      }
    this->Superclass::Z.next_out = buf;
    this->Superclass::Z.avail_out = requestSize;

    do
      {
      if (this->Superclass::Z.avail_in == 0)
        {
        this->Superclass::Z.next_in = this->Superclass::Inbuf;
        this->Superclass::Z.avail_in = static_cast<uInt>(fread(this->Superclass::Inbuf, 1,
            VTK_FOAMFILE_INBUFSIZE, this->Superclass::File));
        if (ferror(this->Superclass::File))
          {
          throw this->StackString() << "Fread failed";
          }
        }
      this->Superclass::ZStatus = inflate(&this->Superclass::Z, Z_NO_FLUSH);
      if (this->Superclass::ZStatus == Z_STREAM_END
#if VTK_FOAMFILE_OMIT_CRCCHECK
      // the dummy CRC function causes data error when finalizing
      // so we have to proceed even when a data error is detected
      || this->Superclass::ZStatus == Z_DATA_ERROR
#endif
      )
        {
        break;
        }
      if (this->Superclass::ZStatus != Z_OK)
        {
        throw this->StackString() << "Inflation failed: "
        << (this->Superclass::Z.msg ? this->Superclass::Z.msg : "");
        }
      } while (this->Superclass::Z.avail_out > 0);

    size = requestSize - this->Superclass::Z.avail_out;
    }
  else
    {
    // not compressed
    size = fread(buf, 1, requestSize, this->Superclass::File);
    }

  if (size <= 0)
    {
    // retain the current location bufPtr_ to the end of the buffer so that
    // getc() returns EOF again when called next time
    return false;
    }
  // size > 0
  // reserve the first byte for getback char
  this->Superclass::BufPtr = this->Superclass::Outbuf + 1;
  this->Superclass::BufEndPtr = this->Superclass::BufPtr + size;
  return true;
}

// get next semantically valid character
int vtkNewFoamFile::NextTokenHead()
{
  for (;;)
    {
    int c;
    while (isspace(c = this->Getc())) // isspace() accepts -1 as EOF
      {
      if (c == '\n')
        {
        ++this->Superclass::LineNumber;
#if VTK_FOAMFILE_RECOGNIZE_LINEHEAD
        this->Superclass::WasNewline = true;
#endif
        }
      }
    if (c == '/')
      {
      if ((c = this->Getc()) == '/')
        {
        while ((c = this->Getc()) != EOF && c != '\n')
          ;
        if (c == EOF)
          {
          return c;
          }
        ++this->Superclass::LineNumber;
#if VTK_FOAMFILE_RECOGNIZE_LINEHEAD
        this->Superclass::WasNewline = true;
#endif
        }
      else if (c == '*')
        {
        for (;;)
          {
          while ((c = this->Getc()) != EOF && c != '*')
            {
            if (c == '\n')
              {
              ++this->Superclass::LineNumber;
              }
            }
          if (c == EOF)
            {
            return c;
            }
          else if ((c = this->Getc()) == '/')
            {
            break;
            }
          this->PutBack(c);
          }
        }
      else
        {
        this->PutBack(c); // may be an EOF
        return '/';
        }
      }
    else // may be an EOF
      {
      return c;
      }
    }
#if defined(__hpux)
  return EOF; // this line should not be executed; workaround for HP-UXia64-aCC
#endif
}

//-----------------------------------------------------------------------------
// class vtkNewFoamIOobject
// holds file handle, file format, name of the object the file holds and
// type of the object.
struct vtkNewFoamIOobject : public vtkNewFoamFile
{
private:
  typedef vtkNewFoamFile Superclass;

public:
  enum fileFormat
    {UNDEFINED, ASCII, BINARY};

private:
  fileFormat Format;
  vtkStdString ObjectName;
  vtkStdString HeaderClassName;
  bool Is13Positions;
  bool IsSinglePrecisionBinary;
  vtkNewFoamError E;

  vtkNewFoamIOobject();
  void ReadHeader(); // defined later
public:
  vtkNewFoamIOobject(const vtkStdString& casePath, const bool isSinglePrecisionBinary) :
    vtkNewFoamFile(casePath), Format(UNDEFINED), Is13Positions(false),
    IsSinglePrecisionBinary(isSinglePrecisionBinary), E()
  {
  }
  ~vtkNewFoamIOobject()
  {
    this->Close();
  }

  bool Open(const vtkStdString& file)
  {
    try
      {
      this->Superclass::Open(file);
      }
    catch(vtkNewFoamError& e)
      {
      this->E = e;
      return false;
      }

    try
      {
      this->ReadHeader();
      }
    catch(vtkNewFoamError& e)
      {
      this->Superclass::Close();
      this->E = e;
      return false;
      }
    return true;
  }

  void Close()
  {
    this->Superclass::Close();
    this->Format = UNDEFINED;
    this->ObjectName.erase();
    this->HeaderClassName.erase();
    this->E.erase();
  }
  fileFormat GetFormat() const
  {
    return this->Format;
  }
  const vtkStdString& GetClassName() const
  {
    return this->HeaderClassName;
  }
  const vtkStdString& GetObjectName() const
  {
    return this->ObjectName;
  }
  const vtkNewFoamError& GetError() const
  {
    return this->E;
  }
  void SetError(const vtkNewFoamError& e)
  {
    this->E = e;
  }
  void SetIs13Positions(const bool is13Positions)
  {
    this->Is13Positions = is13Positions;
  }
  bool GetIs13Positions() const
  {
    return this->Is13Positions;
  }
  bool GetIsSinglePrecisionBinary() const
  {
    return this->IsSinglePrecisionBinary;
  }
};

//-----------------------------------------------------------------------------
// workarounding class for older compilers (gcc-3.3.x and possibly older)
template <typename T> struct vtkNewFoamReadValue
{
public:
  static T ReadValue(vtkNewFoamIOobject &io);
};

VTK_TEMPLATE_SPECIALIZE inline int vtkNewFoamReadValue<int>::ReadValue(vtkNewFoamIOobject& io)
{
  return io.ReadIntValue();
}

VTK_TEMPLATE_SPECIALIZE inline float vtkNewFoamReadValue<float>::ReadValue(vtkNewFoamIOobject& io)
{
  return io.ReadFloatValue();
}

//-----------------------------------------------------------------------------
// class vtkNewFoamEntryValue
// a class that represents a value of a dictionary entry that corresponds to
// its keyword. note that an entry can have more than one value.
struct vtkNewFoamEntryValue : public vtkNewFoamToken
{
public:
  enum uniformTypes
    {UNDEFINED, UNIFORM, NONUNIFORM};
private:
  typedef vtkNewFoamToken Superclass;

  uniformTypes IsUniform;
  bool Managed;
  const vtkNewFoamEntry *UpperEntryPtr;

  vtkNewFoamEntryValue();
  vtkObjectBase *ToVTKObject()
  {
    return this->Superclass::VtkObjectPtr;
  }
  void Clear();
  void ReadList(vtkNewFoamIOobject& io);

public:
  // reads primitive int/float lists
  template <typename listT, typename primitiveT> class listTraits
  {
    listT *Ptr;

  public:
    listTraits() :
      Ptr(listT::New())
    {
    }
    listT *GetPtr()
    {
      return this->Ptr;
    }
    void ReadUniformValues(vtkNewFoamIOobject& io, const int size)
    {
      primitiveT value = vtkNewFoamReadValue<primitiveT>::ReadValue(io);
      for (int i = 0; i < size; i++)
        {
        this->Ptr->SetValue(i, value);
        }
    }
    void ReadAsciiList(vtkNewFoamIOobject& io, const int size)
    {
      for (int i = 0; i < size; i++)
        {
        this->Ptr->SetValue(i, vtkNewFoamReadValue<primitiveT>::ReadValue(io));
        }
    }
    void ReadBinaryList(vtkNewFoamIOobject& io, const int size)
    {
      io.Read(reinterpret_cast<unsigned char *>(this->Ptr->GetPointer(0)), size
          * sizeof(primitiveT));
    }
    void ReadValue(vtkNewFoamIOobject&, vtkNewFoamToken& currToken)
    {
      if (!currToken.Is<primitiveT>())
        {
        throw vtkNewFoamError() << "Expected an integer or a (, found "
        << currToken;
        }
      this->Ptr->InsertNextValue(currToken.To<primitiveT>());
    }
  };

  // reads rank 1 lists of types vector, sphericalTensor, symmTensor
  // and tensor. if isPositions is true it reads Cloud type of data as
  // particle positions. cf. (the positions format)
  // src/lagrangian/basic/particle/particleIO.C
  template <typename listT, typename primitiveT, int nComponents,
      bool isPositions> class vectorListTraits
  {
    listT *Ptr;

  public:
    vectorListTraits() :
      Ptr(listT::New())
    {
      this->Ptr->SetNumberOfComponents(nComponents);
    }
    listT *GetPtr()
    {
      return this->Ptr;
    }
    void ReadUniformValues(vtkNewFoamIOobject& io, const int size)
    {
      io.ReadExpecting('(');
      primitiveT vectorValue[nComponents];
      for (int j = 0; j < nComponents; j++)
        {
        vectorValue[j] = vtkNewFoamReadValue<primitiveT>::ReadValue(io);
        }
      for (int i = 0; i < size; i++)
        {
        this->Ptr->SetTuple(i, vectorValue);
        }
      io.ReadExpecting(')');
      if (isPositions)
        {
        // skip label celli
        vtkNewFoamReadValue<int>::ReadValue(io);
        }
    }
    void ReadAsciiList(vtkNewFoamIOobject& io, const int size)
    {
      for (int i = 0; i < size; i++)
        {
        io.ReadExpecting('(');
        primitiveT *vectorTupleI = this->Ptr->GetPointer(nComponents * i);
        for (int j = 0; j < nComponents; j++)
          {
          vectorTupleI[j] = vtkNewFoamReadValue<primitiveT>::ReadValue(io);
          }
        io.ReadExpecting(')');
        if (isPositions)
          {
          // skip label celli
          vtkNewFoamReadValue<int>::ReadValue(io);
          }
        }
    }
    void ReadBinaryList(vtkNewFoamIOobject& io, const int size)
    {
      if (isPositions) // lagrangian/positions (class Cloud)
        {
        // allocate space along with the larger 1.4 format since the
        // size must be determined at compile-time. The buffer is
        // allocated on the stack in order to avoid leak when an
        // exception is thrown.
        unsigned char buffer[sizeof(double) * (nComponents + 1) + 2
            * sizeof(int)];
        const int realSize = io.GetIsSinglePrecisionBinary()
            ? sizeof(float) : sizeof(double);
        const int nBytes = (io.GetIs13Positions()
            // skip label celli
            ? realSize * nComponents + sizeof(int)
            // skip label celli, label facei and scalar stepFraction
            : realSize * (nComponents + 1) + 2 * sizeof(int));
        for (int i = 0; i < size; i++)
          {
          io.ReadExpecting('(');
          io.Read(buffer, nBytes);
          this->Ptr->SetTuple(i, reinterpret_cast<double *>(buffer));
          io.ReadExpecting(')');
          }
        }
      else
        {
        if (io.GetIsSinglePrecisionBinary())
          {
          io.Read(reinterpret_cast<unsigned char *>(this->Ptr->GetPointer(0)),
              sizeof(float) * nComponents * size);
          }
        else
          {
#if 0
          // reduce calls to io.Read() by using the array space as buffer
          float *destination = this->Ptr->GetPointer(0);
          double *source = reinterpret_cast<double *>(destination);

          int remainingSize = size;
          for (int halfSize = remainingSize / 2; halfSize > 0;
               halfSize = remainingSize / 2)
            {
            const int convertSize = nComponents * halfSize;
            io.Read(reinterpret_cast<unsigned char *>(source),
            sizeof(double) * convertSize);
            for (int i = 0; i < convertSize; i++)
              {
              *destination++ = static_cast<float>(*source++);
              }
            source = reinterpret_cast<double *>(destination);
            remainingSize -= halfSize;
            }
          if (remainingSize > 0)
            {
            // read the last one element that doesn't fit into the
            // remaining space
            double buffer[nComponents];
            io.Read(reinterpret_cast<unsigned char *>(buffer),
                sizeof(double) * nComponents);
            for (int i = 0; i < nComponents; i++)
              {
              *destination++ = static_cast<float>(buffer[i]);
              }
            }
#else
          const int bufferUnit = 32, nDivs = size / bufferUnit,
            unitSize = nComponents * bufferUnit;
          double buffer[nComponents * bufferUnit];

          for (int i = 0; i < nDivs; i++)
            {
            float *destination = this->Ptr->GetPointer(i * unitSize);
            io.Read(reinterpret_cast<unsigned char *>(buffer),
                sizeof(double) * unitSize);
            for (int j = 0; j < unitSize; j++)
              {
              destination[j] = static_cast<float>(buffer[j]);
              }
            }
          const int remainingSize = nComponents * (size % bufferUnit);
          float *destination = this->Ptr->GetPointer(nDivs * unitSize);
          io.Read(reinterpret_cast<unsigned char *>(buffer),
              sizeof(double) * remainingSize);
          for (int j = 0; j < remainingSize; j++)
            {
            destination[j] = static_cast<float>(buffer[j]);
            }
#endif
          }
        }
    }
    void ReadValue(vtkNewFoamIOobject& io, vtkNewFoamToken& currToken)
    {
      if (currToken != '(')
        {
        throw vtkNewFoamError() << "Expected '(', found " << currToken;
        }
      primitiveT v[nComponents];
      for (int j = 0; j < nComponents; j++)
        {
        v[j] = vtkNewFoamReadValue<primitiveT>::ReadValue(io);
        }
      this->Ptr->InsertNextTuple(v);
      io.ReadExpecting(')');
    }
  };

  vtkNewFoamEntryValue(const vtkNewFoamEntry *upperEntryPtr) :
    vtkNewFoamToken(), IsUniform(UNDEFINED), Managed(true),
        UpperEntryPtr(upperEntryPtr)
  {
  }
  vtkNewFoamEntryValue(vtkNewFoamEntryValue&, const vtkNewFoamEntry *);
  ~vtkNewFoamEntryValue()
  {
    this->Clear();
  }

  void SetEmptyList()
  {
    this->Clear();
    this->IsUniform = NONUNIFORM;
    this->Superclass::Type = EMPTYLIST;
  }
  uniformTypes GetIsUniform() const
  {
    return this->IsUniform;
  }
  void SetIsUniform(const uniformTypes isUniform)
  {
    this->IsUniform = isUniform;
  }
  void Read(vtkNewFoamIOobject& io);
  void ReadDictionary(vtkNewFoamIOobject& io, const vtkNewFoamToken& firstKeyword);
  const vtkIntArray& LabelList() const
  {
    return *this->Superclass::LabelListPtr;
  }
  vtkIntArray& LabelList()
  {
    return *this->Superclass::LabelListPtr;
  }
  const vtkNewFoamIntVectorVector& LabelListList() const
  {
    return *this->Superclass::LabelListListPtr;
  }
  const vtkFloatArray& ScalarList() const
  {
    return *this->Superclass::ScalarListPtr;
  }
  vtkFloatArray& ScalarList()
  {
    return *this->Superclass::ScalarListPtr;
  }
  const vtkFloatArray& VectorList() const
  {
    return *this->Superclass::VectorListPtr;
  }
  const vtkNewFoamDict& Dictionary() const
  {
    return *this->Superclass::DictPtr;
  }
  vtkNewFoamDict& Dictionary()
  {
    return *this->Superclass::DictPtr;
  }

  void *Ptr()
  {
    this->Managed = false; // returned pointer will not be deleted by the d'tor
    // all list pointers are in a single union
    return (void *)this->Superclass::LabelListPtr;
  }

  vtkStdString ToStdString() const
  {
    return this->Superclass::IsWordOrString()
        || this->Superclass::GetType() == vtkNewFoamToken::IDENTIFIER
        ? this->Superclass::ToStdString() : vtkStdString();
  }
  vtkStdString ToWord() const
  {
    return this->Superclass::Type == WORD ? this->Superclass::ToStdString()
        : vtkStdString();
  }
  vtkStdString ToString() const
  {
    return this->Superclass::Type == STRING ? this->Superclass::ToStdString()
        : vtkStdString();
  }
  float ToFloat() const
  {
    return this->Superclass::Type == SCALAR || this->Superclass::Type == LABEL ? this->Superclass::To<float>()
        : 0.0F;
  }
  double ToDouble() const
  {
    return this->Superclass::Type == SCALAR || this->Superclass::Type == LABEL ? this->Superclass::To<double>()
        : 0.0;
  }
  int ToInt() const
  {
    return this->Superclass::Type == LABEL ? this->Superclass::To<int>() : 0;
  }

  // the following two are for an exceptional expression of
  // `LABEL{LABELorSCALAR}' without type prefix (e. g. `2{-0}' in
  // mixedRhoE B.C. in rhopSonicFoam/shockTube)
  void MakeLabelList(const int labelValue, const int size)
  {
    this->Superclass::LabelListPtr = vtkIntArray::New();
    this->Superclass::Type = LABELLIST;
    this->Superclass::LabelListPtr->SetNumberOfValues(size);
    for (int i = 0; i < size; i++)
      {
      this->Superclass::LabelListPtr->SetValue(i, labelValue);
      }
  }
  void MakeScalarList(const float scalarValue, const int size)
  {
    this->Superclass::ScalarListPtr = vtkFloatArray::New();
    this->Superclass::Type = SCALARLIST;
    this->Superclass::ScalarListPtr->SetNumberOfValues(size);
    for (int i = 0; i < size; i++)
      {
      this->Superclass::ScalarListPtr->SetValue(i, scalarValue);
      }
  }

  // reads dimensionSet
  void ReadDimensionSet(vtkNewFoamIOobject& io)
  {
    const int nDims = 7;
    this->Superclass::LabelListPtr = vtkIntArray::New();
    this->Superclass::Type = LABELLIST;
    this->Superclass::LabelListPtr->SetNumberOfValues(nDims);
    for (int dimI = 0; dimI < nDims; dimI++)
      {
      this->Superclass::LabelListPtr->SetValue(dimI, vtkNewFoamReadValue<int>::ReadValue(io));
      }
    io.ReadExpecting(']');
  }

  template <vtkNewFoamToken::tokenType listType, typename traitsT> void ReadNonuniformList(
      vtkNewFoamIOobject& io);

  // reads a list of labelLists. requires size prefix of the listList
  // to be present. size of each sublist must also be present in the
  // stream if the format is binary.
  void ReadLabelListList(vtkNewFoamIOobject& io)
  {
    vtkNewFoamToken currToken;
    if (!io.Read(currToken))
      {
      throw vtkNewFoamError() << "Unexpected EOF";
      }
    if (currToken.GetType() == vtkNewFoamToken::LABEL)
      {
      const int sizeI = currToken.To<int>();
      if (sizeI < 0)
        {
        throw vtkNewFoamError() << "List size must not be negative: size = "
        << sizeI;
        }
      // gives initial guess for list size
      this->Superclass::LabelListListPtr = new vtkNewFoamIntVectorVector(sizeI, 4 * sizeI);
      this->Superclass::Type = LABELLISTLIST;
      io.ReadExpecting('(');
      int bodyI = 0;
      for (int i = 0; i < sizeI; i++)
        {
        if (!io.Read(currToken))
          {
          throw vtkNewFoamError() << "Unexpected EOF";
          }
        if (currToken.GetType() == vtkNewFoamToken::LABEL)
          {
          const int sizeJ = currToken.To<int>();
          if (sizeJ < 0)
            {
            throw vtkNewFoamError() << "List size must not be negative: size = "
            << sizeJ;
            }
          if (bodyI + sizeJ > this->Superclass::LabelListListPtr->GetBodySize())
            {
            const int newSize =
                this->Superclass::LabelListListPtr->GetBodySize() + sizeJ;
            this->Superclass::LabelListListPtr->ResizeBody(newSize);
            }
          int *listI = this->Superclass::LabelListListPtr->SetIndex(i, bodyI);
          if (io.GetFormat() == vtkNewFoamIOobject::ASCII)
            {
            io.ReadExpecting('(');
            for (int j = 0; j < sizeJ; j++)
              {
              listI[j] = vtkNewFoamReadValue<int>::ReadValue(io);
              }
            io.ReadExpecting(')');
            }
          else
            {
            if (sizeJ > 0) // avoid invalid reference to labelListI.at(0)
              {
              io.ReadExpecting('(');
              io.Read(reinterpret_cast<unsigned char*>(listI), sizeJ
                  * sizeof(int));
              io.ReadExpecting(')');
              }
            }
          bodyI += sizeJ;
          }
        else if (currToken == '(')
          {
          this->Superclass::LabelListListPtr->SetIndex(i, bodyI);
          while (io.Read(currToken) && currToken != ')')
            {
            if (currToken.GetType() != vtkNewFoamToken::LABEL)
              {
              throw vtkNewFoamError() << "Expected an integer, found "
              << currToken;
              }
            if (bodyI >= this->LabelListListPtr->GetBodySize())
              {
              const int newSize =
                  this->Superclass::LabelListListPtr->GetBodySize() + 1;
              this->Superclass::LabelListListPtr->ResizeBody(newSize);
              }
            this->Superclass::LabelListListPtr
            ->SetValue(bodyI++, currToken.To<int>());
            }
          }
        else
          {
          throw vtkNewFoamError() << "Expected integer or '(', found "
          << currToken;
          }
        }
      // set the next index of the last element to calculate the last
      // subarray size
      this->Superclass::LabelListListPtr->SetIndex(sizeI, bodyI);
      // shrink to the actually used size
      this->Superclass::LabelListListPtr->ResizeBody(bodyI);
      io.ReadExpecting(')');
      }
    else
      {
      throw vtkNewFoamError() << "Expected integer, found " << currToken;
      }
  }

  bool ReadField(vtkNewFoamIOobject& io)
  {
    try
      {
      // lagrangian labels (cf. gnemdFoam/nanoNozzle)
      if(io.GetClassName() == "labelField")
        {
        this->ReadNonuniformList<LABELLIST, listTraits<vtkIntArray, int> >(io);
        }
      // lagrangian scalars

      else if(io.GetClassName() == "scalarField")
        {
        this->ReadNonuniformList<SCALARLIST, listTraits<vtkFloatArray, float> >(
            io);
        }
      else if(io.GetClassName() == "sphericalTensorField")
        {
        this->ReadNonuniformList<VECTORLIST,
        vectorListTraits<vtkFloatArray, float, 1, false> >(io);
        }
      // polyMesh/points, lagrangian vectors

      else if(io.GetClassName() == "vectorField")
        {
        this->ReadNonuniformList<VECTORLIST,
        vectorListTraits<vtkFloatArray, float, 3, false> >(io);
        }
      else if(io.GetClassName() == "symmTensorField")
        {
        this->ReadNonuniformList<VECTORLIST,
        vectorListTraits<vtkFloatArray, float, 6, false> >(io);
        }
      else if(io.GetClassName() == "tensorField")
        {
        this->ReadNonuniformList<VECTORLIST,
        vectorListTraits<vtkFloatArray, float, 9, false> >(io);
        }
      else
        {
        throw vtkNewFoamError() << "Non-supported field type "
        << io.GetClassName();
        }
      }
    catch(vtkNewFoamError& e)
      {
      io.SetError(e);
      return false;
      }
    return true;
  }
};

// specialization for reading double precision binary into vtkFloatArray.
// Must precede ReadNonuniformList() below (HP-UXia64-aCC).
VTK_TEMPLATE_SPECIALIZE
void vtkNewFoamEntryValue::listTraits<vtkFloatArray, float>::ReadBinaryList(
    vtkNewFoamIOobject& io, const int size)
{
  if(io.GetIsSinglePrecisionBinary())
    {
    io.Read(reinterpret_cast<unsigned char *>(this->Ptr->GetPointer(0)),
        sizeof(float) * size);
    }
  else
    {
#if 0
    // reduce calls to io.Read() by using the array space as buffer
    float *destination = this->Ptr->GetPointer(0);
    double *source = reinterpret_cast<double *>(destination);

    int remainingSize = size;
    for (int halfSize = size / 2; halfSize > 0; halfSize = remainingSize / 2)
      {
      io.Read(reinterpret_cast<unsigned char *>(source),
          sizeof(double) * halfSize);
      for (int i = 0; i < halfSize; i++)
        {
        *destination++ = static_cast<float>(*source++);
        }
      source = reinterpret_cast<double *>(destination);
      remainingSize -= halfSize;
      }
    if (remainingSize > 0)
      {
      // read the last one element that doesn't fit into the remaining space
      double buffer;
      io.Read(reinterpret_cast<unsigned char *>(&buffer), sizeof(double));
      *destination = static_cast<float>(buffer);
      }
#else
    const int bufferUnit = 32, nDivs = size / bufferUnit;
    double buffer[bufferUnit];

    for (int i = 0; i < nDivs; i++)
      {
      float *destination = this->Ptr->GetPointer(i * bufferUnit);
      io.Read(reinterpret_cast<unsigned char *>(buffer),
          sizeof(double) * bufferUnit);
      for (int j = 0; j < bufferUnit; j++)
        {
        destination[j] = static_cast<float>(buffer[j]);
        }
      }
    const int remainingSize = size % bufferUnit;
    float *destination = this->Ptr->GetPointer(nDivs * bufferUnit);
    io.Read(reinterpret_cast<unsigned char *>(buffer),
        sizeof(double) * remainingSize);
    for (int j = 0; j < remainingSize; j++)
      {
      destination[j] = static_cast<float>(buffer[j]);
      }
#endif
    }
}

// generic reader for nonuniform lists. requires size prefix of the
// list to be present in the stream if the format is binary.
template <vtkNewFoamToken::tokenType listType, typename traitsT>
void vtkNewFoamEntryValue::ReadNonuniformList(vtkNewFoamIOobject& io)
{
  vtkNewFoamToken currToken;
  if (!io.Read(currToken))
    {
    throw vtkNewFoamError() << "Unexpected EOF";
    }
  traitsT list;
  this->Superclass::Type = listType;
  this->Superclass::VtkObjectPtr = list.GetPtr();
  if (currToken.Is<int>())
    {
    const int size = currToken.To<int>();
    if (size < 0)
      {
      throw vtkNewFoamError() << "List size must not be negative: size = " << size;
      }
    list.GetPtr()->SetNumberOfTuples(size);
    if (io.GetFormat() == vtkNewFoamIOobject::ASCII)
      {
      if (!io.Read(currToken))
        {
        throw vtkNewFoamError() << "Unexpected EOF";
        }
      // some objects have lists with only one element enclosed by {}
      // e. g. simpleFoam/pitzDaily3Blocks/constant/polyMesh/faceZones
      if (currToken == '{')
        {
        list.ReadUniformValues(io, size);
        io.ReadExpecting('}');
        return;
        }
      else if (currToken != '(')
        {
        throw vtkNewFoamError() << "Expected '(', found " << currToken;
        }
      list.ReadAsciiList(io, size);
      io.ReadExpecting(')');
      }
    else
      {
      if (size > 0)
        {
        // read parentheses only when size > 0
        io.ReadExpecting('(');
        list.ReadBinaryList(io, size);
        io.ReadExpecting(')');
        }
      }
    }
  else if (currToken == '(')
    {
    while (io.Read(currToken) && currToken != ')')
      {
      list.ReadValue(io, currToken);
      }
    list.GetPtr()->Squeeze();
    }
  else
    {
    throw vtkNewFoamError() << "Expected integer or '(', found " << currToken;
    }
}

//-----------------------------------------------------------------------------
// class vtkNewFoamKeyword
// a class that handles regex
#ifdef VTK_FOAMFILE_HAVE_REGEX
// use the system POSIX regex if it was found
struct vtkNewFoamKeyword : public vtkStdString
{
private:
  typedef vtkStdString Superclass;
  regex_t *Preg;
  void operator=(const vtkNewFoamKeyword &); // not implemented

  void Clear()
  {
    if (this->Preg != NULL)
      {
      regfree(this->Preg);
      delete this->Preg;
      this->Preg = NULL;
      }
  }

  int Compile(const vtkStdString &keyword)
  {
    this->Preg = new regex_t;
    return regcomp(this->Preg, keyword.c_str(), REG_EXTENDED);
  }

public:
  vtkNewFoamKeyword()
    : vtkStdString(), Preg(NULL)
  {
  }
  vtkNewFoamKeyword(const vtkNewFoamKeyword &keyword)
    : vtkStdString(keyword), Preg(NULL)
  {
    if (keyword.Preg != NULL)
      {
      // does not handle error assuming the compilation always succeed
      // since exception should not be thrown from within the constructor
      this->Compile(keyword);
      }
  }
  ~vtkNewFoamKeyword()
  {
    this->Clear();
  }

  void SetKeyword(const vtkNewFoamToken &keywordToken, const bool isRegEx)
  {
    this->Clear();

    this->Superclass::operator=(keywordToken.ToStdString());
    if (keywordToken.GetType() == vtkNewFoamToken::STRING && isRegEx)
      {
      const int ret = this->Compile(keywordToken.ToStdString());
      if (ret != 0)
        {
        const int msgSize = regerror(ret, this->Preg, NULL, 0);
        char *msgBuf = new char[msgSize + 1];
        regerror(ret, this->Preg, msgBuf, msgSize + 1);
        vtkStdString errString = msgBuf;
        delete [] msgBuf;
        this->Clear();
        throw vtkNewFoamError() << "regular expression error: " << errString;
        }
      }
  }

  bool RegExMatch(const vtkStdString &string) const
  {
    if (this->Preg == NULL)
      {
      return false;
      }

    regmatch_t pmatch;
    const int ret = regexec(this->Preg, string.c_str(), 1, &pmatch, 0);
    return ret == 0 && pmatch.rm_so == 0
        && static_cast<size_t>(pmatch.rm_eo) == string.length();
  }
};
#elif defined(_MSC_VER) && (_MSC_FULL_VER > 150030729 \
    || _MSC_FULL_VER == 150030729 && _MSC_BUILD >= 1)
// if VS2008SP1 or above, use std::tr1::regex
#include <regex>
struct vtkNewFoamKeyword : public vtkStdString
{
private:
  typedef vtkStdString Superclass;
  vtkstd::tr1::regex *Regex;
  void operator=(const vtkNewFoamKeyword &); // not implemented

  void Clear()
  {
    delete this->Regex;
    this->Regex = NULL;
  }

public:
  vtkNewFoamKeyword()
    : vtkStdString(), Regex(NULL)
  {
  }
  vtkNewFoamKeyword(const vtkNewFoamKeyword &keyword)
    : vtkStdString(keyword), Regex(NULL)
  {
    if (keyword.Regex != NULL)
      {
      this->Regex = new vtkstd::tr1::regex(*keyword.Regex);
      }
  }
  ~vtkNewFoamKeyword()
  {
    this->Clear();
  }

  void SetKeyword(const vtkNewFoamToken &keywordToken, const bool isRegEx)
  {
    this->Clear();

    this->Superclass::operator=(keywordToken.ToStdString());
    if (keywordToken.GetType() == vtkNewFoamToken::STRING && isRegEx)
      {
      try
        {
        this->Regex = new vtkstd::tr1::regex(
            keywordToken.ToStdString(), vtkstd::tr1::regex::extended);
        }
      catch (const vtkstd::tr1::regex_error)
        {
        this->Clear();
        throw vtkNewFoamError() << "regular expression error";
        }
      }
  }

  bool RegExMatch(const vtkStdString &string) const
  {
    if (this->Regex == NULL)
      {
      return false;
      }

    return regex_match(string, *Regex);
  }
};
#else
// if none is found, use vtksys::RegularExpression as the fallback option
struct vtkNewFoamKeyword : public vtkStdString
{
private:
  typedef vtkStdString Superclass;
  vtksys::RegularExpression *Regex;
  void operator=(const vtkNewFoamKeyword &); // not implemented

  void Clear()
  {
    delete this->Regex;
    this->Regex = NULL;
  }

public:
  vtkNewFoamKeyword()
    : vtkStdString(), Regex(NULL)
  {
  }
  vtkNewFoamKeyword(const vtkNewFoamKeyword &keyword)
    : vtkStdString(keyword), Regex(NULL)
  {
    if (keyword.Regex != NULL)
      {
      this->Regex = new vtksys::RegularExpression(*keyword.Regex);
      }
  }
  ~vtkNewFoamKeyword()
  {
    this->Clear();
  }

  void SetKeyword(const vtkNewFoamToken &keywordToken, const bool isRegEx)
  {
    this->Clear();

    this->Superclass::operator=(keywordToken.ToStdString());
    if (keywordToken.GetType() == vtkNewFoamToken::STRING && isRegEx)
      {
      this->Regex = new vtksys::RegularExpression(
          keywordToken.ToStdString().c_str());
      if (!this->Regex->is_valid())
        {
        this->Clear();
        throw vtkNewFoamError() << "regular expression error";
        }
      }
  }

  bool RegExMatch(const vtkStdString &string) const
  {
    if (this->Regex == NULL)
      {
      return false;
      }

    return this->Regex->find(string) && this->Regex->start() == 0
        && this->Regex->end() == string.length();
  }
};
#endif

//-----------------------------------------------------------------------------
// class vtkNewFoamEntry
// a class that represents an entry of a dictionary. note that an
// entry can have more than one value.
struct vtkNewFoamEntry : public vtkstd::vector<vtkNewFoamEntryValue*>
{
private:
  typedef vtkstd::vector<vtkNewFoamEntryValue*> Superclass;
  vtkNewFoamKeyword Keyword;
  vtkNewFoamDict *UpperDictPtr;

  vtkNewFoamEntry();

public:
  vtkNewFoamEntry(vtkNewFoamDict *upperDictPtr) :
    Keyword(), UpperDictPtr(upperDictPtr)
  {
  }
  vtkNewFoamEntry(const vtkNewFoamEntry& entry, vtkNewFoamDict *upperDictPtr) :
    Superclass(entry.size()), Keyword(entry.GetKeyword()),
        UpperDictPtr(upperDictPtr)
  {
    for (size_t valueI = 0; valueI < entry.size(); valueI++)
      {
      this->Superclass::operator[](valueI) = new vtkNewFoamEntryValue(*entry[valueI], this);
      }
  }

  ~vtkNewFoamEntry()
  {
    this->Clear();
  }

  void Clear()
  {
    for (size_t i = 0; i < this->Superclass::size(); i++)
      {
      delete this->Superclass::operator[](i);
      }
    this->Superclass::clear();
  }
  const vtkNewFoamKeyword& GetKeyword() const
  {
    return this->Keyword;
  }
  void SetKeyword(const vtkNewFoamToken& keyword, const bool isRegExKeyword)
  {
    this->Keyword.SetKeyword(keyword, isRegExKeyword);
  }
  const vtkNewFoamEntryValue& FirstValue() const
  {
    return *this->Superclass::operator[](0);
  }
  vtkNewFoamEntryValue& FirstValue()
  {
    return *this->Superclass::operator[](0);
  }
  const vtkIntArray& LabelList() const
  {
    return this->FirstValue().LabelList();
  }
  vtkIntArray& LabelList()
  {
    return this->FirstValue().LabelList();
  }
  const vtkNewFoamIntVectorVector& LabelListList() const
  {
    return this->FirstValue().LabelListList();
  }
  const vtkFloatArray& ScalarList() const
  {
    return this->FirstValue().ScalarList();
  }
  vtkFloatArray& ScalarList()
  {
    return this->FirstValue().ScalarList();
  }
  const vtkFloatArray& VectorList() const
  {
    return this->FirstValue().VectorList();
  }
  const vtkNewFoamDict& Dictionary() const
  {
    return this->FirstValue().Dictionary();
  }
  vtkNewFoamDict& Dictionary()
  {
    return this->FirstValue().Dictionary();
  }
  void *Ptr()
  {
    return this->FirstValue().Ptr();
  }
  const vtkNewFoamDict *GetUpperDictPtr() const
  {
    return this->UpperDictPtr;
  }

  vtkStdString ToStdString() const
  {
    return this->Superclass::size() > 0 ? this->FirstValue().ToStdString() : vtkStdString();
  }
  vtkStdString ToWord() const
  {
    return this->Superclass::size() > 0 ? this->FirstValue().ToWord() : vtkStdString();
  }
  vtkStdString ToString() const
  {
    return this->Superclass::size() > 0 ? this->FirstValue().ToString() : vtkStdString();
  }
  float ToFloat() const
  {
    return this->Superclass::size() > 0 ? this->FirstValue().ToFloat() : 0.0F;
  }
  double ToDouble() const
  {
    return this->Superclass::size() > 0 ? this->FirstValue().ToDouble() : 0.0;
  }
  int ToInt() const
  {
    return this->Superclass::size() > 0 ? this->FirstValue().ToInt() : 0;
  }

  void ReadDictionary(vtkNewFoamIOobject& io)
  {
    this->Superclass::push_back(new vtkNewFoamEntryValue(this));
    this->Superclass::back()->ReadDictionary(io, vtkNewFoamToken());
  }

  // read values of an entry
  void Read(vtkNewFoamIOobject& io);
};

//-----------------------------------------------------------------------------
// class vtkNewFoamDict
// a class that holds a FoamFile data structure
struct vtkNewFoamDict : public vtkstd::vector<vtkNewFoamEntry*>
{
private:
  typedef vtkstd::vector<vtkNewFoamEntry*> Superclass;

  vtkNewFoamToken Token;
  const vtkNewFoamDict *UpperDictPtr;

  vtkNewFoamDict(const vtkNewFoamDict &);

public:
  vtkNewFoamDict(const vtkNewFoamDict *upperDictPtr = NULL) :
    Superclass(), Token(), UpperDictPtr(upperDictPtr)
  {
  }
  vtkNewFoamDict(const vtkNewFoamDict& dict, const vtkNewFoamDict *upperDictPtr) :
    Superclass(dict.size()), Token(), UpperDictPtr(upperDictPtr)
  {
    if (dict.GetType() == vtkNewFoamToken::DICTIONARY)
      {
      for (size_t entryI = 0; entryI < dict.size(); entryI++)
        {
        this->operator[](entryI) = new vtkNewFoamEntry(*dict[entryI], this);
        }
      }
  }

  ~vtkNewFoamDict()
  {
    if (this->Token.GetType() == vtkNewFoamToken::UNDEFINED)
      {
      for (size_t i = 0; i < this->Superclass::size(); i++)
        {
        delete this->operator[](i);
        }
      }
  }

  vtkNewFoamToken::tokenType GetType() const
  {
    return this->Token.GetType() == vtkNewFoamToken::UNDEFINED ? vtkNewFoamToken::DICTIONARY
        : this->Token.GetType();
  }
  const vtkNewFoamToken &GetToken() const
  {
    return this->Token;
  }
  const vtkNewFoamDict *GetUpperDictPtr() const
  {
    return this->UpperDictPtr;
  }

  // lookup by literal search
  vtkNewFoamEntry *LookupLiteral(const vtkStdString& keyword) const
  {
    if (this->Token.GetType() == vtkNewFoamToken::UNDEFINED)
      {
      for (size_t entryI = 0; entryI < this->Superclass::size(); entryI++)
        {
        // attempt literal match
        if (this->operator[](entryI)->GetKeyword() == keyword)
          {
          // found
          return this->operator[](entryI);
          }
        }
      }

    // not found
    return NULL;
  }

  // lookup by literal and regex search
  vtkNewFoamEntry *LookupRegEx(const vtkStdString& keyword) const
  {
    if (this->Token.GetType() == vtkNewFoamToken::UNDEFINED)
      {
      // first attempt a literal match regardless of the keywords being
      // matched containing regular expressions or not
      vtkNewFoamEntry *entry = this->LookupLiteral(keyword);
      if (entry == NULL)
        {
        // regex is matched in descending order
        for (int entryI = static_cast<int>(this->Superclass::size()) - 1;
            entryI >= 0; entryI--)
          {
          // attempt regex match
          if (this->operator[](entryI)->GetKeyword().RegExMatch(keyword))
            {
            // found
            entry = this->operator[](entryI);
            break;
            }
          }
        }
      return entry;
      }

    return NULL;
  }

  // return NULL if the found entry contains no value
  vtkNewFoamEntry *Lookup(const vtkStdString &keyword) const
  {
    vtkNewFoamEntry *entry = this->LookupRegEx(keyword);
    return (entry != NULL && entry->size() > 0) ? entry : NULL;
  }

  // reads a FoamFile or a subdictionary. if the stream to be read is
  // a subdictionary the preceding '{' is assumed to have already been
  // thrown away.
  bool Read(vtkNewFoamIOobject& io, const bool isSubDictionary = false,
      const vtkNewFoamToken& firstToken = vtkNewFoamToken())
  {
    try
      {
      bool isRegExKeyword = true;
      vtkNewFoamToken currToken;
      if (firstToken.GetType() == vtkNewFoamToken::UNDEFINED)
        {
        // read the first token
        if (!io.Read(currToken))
          {
          throw vtkNewFoamError() << "Unexpected EOF";
          }

        if (isSubDictionary)
          {
          // the following if clause is for an exceptional expression
          // of `LABEL{LABELorSCALAR}' without type prefix
          // (e. g. `2{-0}' in mixedRhoE B.C. in
          // rhopSonicFoam/shockTube)
          if (currToken.GetType() == vtkNewFoamToken::LABEL
              || currToken.GetType() == vtkNewFoamToken::SCALAR)
            {
            this->Token = currToken;
            io.ReadExpecting('}');
            return true;
            }
          // return as empty dictionary
          else if (currToken == '}')
            {
            return true;
            }
          }
        else
          {
          // list of dictionaries is read as a usual dictionary with regex
          // turned off: polyMesh/boundary, point/face/cell-Zones
          if (currToken.GetType() == vtkNewFoamToken::LABEL)
            {
            io.ReadExpecting('(');
            if (currToken.To<int>() > 0)
              {
              if (!io.Read(currToken))
                {
                throw vtkNewFoamError() << "Unexpected EOF";
                }
              // continue to read as a usual dictionary with regex turned off
              isRegExKeyword = false;
              }
            else // return as empty dictionary
              {
              io.ReadExpecting(')');
              return true;
              }
            }
          // some boundary files does not have the number of boundary
          // patches (e.g. settlingFoam/tank3D). In this case the file is
          // needed to be explicitly read as a dictionary.
          else if (currToken == '('
              && io.GetClassName() == "polyBoundaryMesh") // polyMesh/boundary
            {
            if (!io.Read(currToken)) // read the first keyword
              {
              throw vtkNewFoamError() << "Unexpected EOF";
              }
            if (currToken == ')') // return as empty dictionary
              {
              return true;
              }
            // continue to read as a usual dictionary with regex turned off
            isRegExKeyword = false;
            }
          }
        }
      // if firstToken is given as either a word or a string read the
      // following stream as subdictionary to constitute a list of
      // dictionaries. Both word and string are allowed as keyword,
      // but string does not work as regex.
      else if (firstToken.IsWordOrString())
        {
        this->Superclass::push_back(new vtkNewFoamEntry(this));
        this->Superclass::back()->SetKeyword(firstToken, isRegExKeyword);
        this->Superclass::back()->ReadDictionary(io);
        if (!io.Read(currToken) || currToken == '}' || currToken == ')')
          {
          return true;
          }
        }
      else // quite likely an identifier
        {
        currToken = firstToken;
        }

      if (currToken == ';' || currToken.IsWordOrString()
          || currToken.GetType() == vtkNewFoamToken::IDENTIFIER)
        {
        // general dictionary
        do
          {
          if (currToken.IsWordOrString())
            {
            vtkNewFoamEntry *previousEntry;
            if (isRegExKeyword && (previousEntry
                = this->LookupLiteral(currToken.ToStdString())) != NULL)
              {
              if (io.GetInputMode() == vtkNewFoamFile::INPUT_MODE_MERGE)
                {
                if (previousEntry->FirstValue().GetType()
                    == vtkNewFoamToken::DICTIONARY)
                  {
                  io.ReadExpecting('{');
                  previousEntry->FirstValue().Dictionary().Read(io, true);
                  }
                else
                  {
                  previousEntry->Clear();
                  previousEntry->Read(io);
                  }
                }
              else if (io.GetInputMode() == vtkNewFoamFile::INPUT_MODE_OVERWRITE)
                {
                previousEntry->Clear();
                previousEntry->Read(io);
                }
              else if (io.GetInputMode() == vtkNewFoamFile::INPUT_MODE_PROTECT)
                {
                // the contents of the entry is discarded
                vtkNewFoamEntry entry(this);
                entry.Read(io);
                }
              else // INPUT_MODE_ERROR || INPUT_MODE_WARN
                {
                // "#inputMode warn" just doesn't skip the duplicated entry but
                // in fact stops reading the remaining dictionary
                throw vtkNewFoamError() << "Found duplicated entries with keyword "
                  << currToken.ToStdString();
                }
              }
            else
              {
              this->Superclass::push_back(new vtkNewFoamEntry(this));
              this->Superclass::back()->SetKeyword(currToken, isRegExKeyword);
              this->Superclass::back()->Read(io);
              }

            // the "FoamFile" keyword of type word must be searched by exact
            // match (no need for regex search)
            if(currToken == "FoamFile")
              {
              // delete the FoamFile header subdictionary entry
              delete this->Superclass::back();
              this->Superclass::pop_back();
              }
            // the "include" keyword of type word or string must be searched by
            // exact match (no need for regex search)
            else if (currToken.ToStdString() == "include")
              {
              // include the named file. The name must be of type string.
              // Exiting the included file at EOF will be handled automatically
              // by vtkNewFoamFile::closeIncludedFile()
                if (this->Superclass::back()->size() == 0
                    || this->Superclass::back()->FirstValue().GetType()
                    != vtkNewFoamToken::STRING)
                {
                throw vtkNewFoamError()
                << "Expected string as the file name to be included, found "
                << this->Superclass::back()->size() == 0 ? vtkNewFoamToken()
                : this->Superclass::back()->FirstValue();
                }
              const vtkStdString includeFileName(
                  this->Superclass::back()->ToString());
              delete this->Superclass::back();
              this->Superclass::pop_back();
              io.IncludeFile(includeFileName, io.GetFilePath(), false);
              }
            }
          else if (currToken.GetType() == vtkNewFoamToken::IDENTIFIER)
            {
            // substitute identifier
            const vtkStdString identifier(currToken.ToStdString());

            for (const vtkNewFoamDict *uDictPtr = this; uDictPtr != NULL;)
              {
              const vtkNewFoamEntry *identifiedEntry
                  = uDictPtr->LookupLiteral(identifier); // do a regex search
              if (identifiedEntry != NULL)
                {
                if (identifiedEntry->size() == 0
                    || identifiedEntry->FirstValue().GetType()
                    != vtkNewFoamToken::DICTIONARY)
                  {
                  throw vtkNewFoamError()
                  << "Expected dictionary for substituting entry "
                  << identifier;
                  }
                const vtkNewFoamDict& identifiedDict
                = identifiedEntry->FirstValue().Dictionary();
                for (size_t entryI = 0; entryI < identifiedDict.size(); entryI++)
                  {
                  // I think #inputMode handling should be done here
                  // as well, but the genuine FoamFile parser for OF
                  // 1.5 does not seem to be doing it.
                  this->Superclass::push_back(
                      new vtkNewFoamEntry(*identifiedDict[entryI], this));
                  }
                break;
                }
              else
                {
                uDictPtr = uDictPtr->GetUpperDictPtr();
                if (uDictPtr == NULL)
                  {
                  // if no entry with the identifier as keyowrd is found
                  // interpret the identifier as keyword
                  this->Superclass::push_back(new vtkNewFoamEntry(this));
                  vtkNewFoamToken keywordToken;
                  keywordToken.SetWord("$" + identifier);
                  this->Superclass::back()->SetKeyword(keywordToken, isRegExKeyword);
                  this->Superclass::back()->Read(io);
                  }
                }
              }
            }
          // skip empty entry only with ';'
          }while (io.Read(currToken) && (currToken.IsWordOrString()
                || currToken.GetType() == vtkNewFoamToken::IDENTIFIER
                || currToken == ';'));

        if (currToken.GetType() == vtkNewFoamToken::TOKEN_ERROR || currToken == '}'
            || currToken == ')')
          {
          return true;
          }
        throw vtkNewFoamError()
        << "Expected keyword, closing brace, ';' or EOF, found " << currToken;
        }
      throw vtkNewFoamError() << "Expected keyword or identifier, found "
      << currToken;
      }
    catch (vtkNewFoamError& e)
      {
      if (isSubDictionary)
        {
        throw;
        }
      else
        {
        io.SetError(e);
        return false;
        }
      }
  }
};

void vtkNewFoamIOobject::ReadHeader()
{
  // the "FoamFile" keyword must be searched by exact match (no need
  // for regex search)
  this->Superclass::ReadExpecting("FoamFile");
  this->Superclass::ReadExpecting('{');

  vtkNewFoamDict headerDict;
  // throw exception in case of error
  headerDict.Read(*this, true, vtkNewFoamToken());

  const vtkNewFoamEntry *formatEntry = headerDict.Lookup("format");
  if (formatEntry == NULL || !formatEntry->FirstValue().IsWordOrString())
    {
    throw vtkNewFoamError()
    << "valid format entry (binary/ascii) not found in FoamFile header";
    }
  // case does matter (e. g. "BINARY" is treated as ascii)
  // cf. src/OpenFOAM/db/IOstreams/IOstreams/IOstream.C
  this->Format = (formatEntry->ToStdString() == "binary" ? BINARY : ASCII);

  const vtkNewFoamEntry *classEntry = headerDict.Lookup("class");
  if (classEntry == NULL || !classEntry->FirstValue().IsWordOrString())
    {
    throw vtkNewFoamError() << "valid class name not found in FoamFile header";
    }
  this->HeaderClassName = classEntry->ToStdString();

  const vtkNewFoamEntry *objectEntry = headerDict.Lookup("object");
  if (objectEntry == NULL || !objectEntry->FirstValue().IsWordOrString())
    {
    throw vtkNewFoamError() << "valid object name not found in FoamFile header";
    }
  this->ObjectName = objectEntry->ToStdString();
}

vtkNewFoamEntryValue::vtkNewFoamEntryValue(
    vtkNewFoamEntryValue& value, const vtkNewFoamEntry *upperEntryPtr) :
  vtkNewFoamToken(value), IsUniform(value.GetIsUniform()), Managed(true),
      UpperEntryPtr(upperEntryPtr)
{
  switch (this->Superclass::Type)
    {
    case VECTORLIST:
#if vtksys_DATE_STAMP_FULL >= 20080620
        {
        vtkFloatArray *fa = vtkFloatArray::SafeDownCast(value.ToVTKObject());
        if(fa->GetNumberOfComponents() == 6)
          {
          // create deepcopies for vtkObjects to avoid duplicated
          // mainpulation of symmTensor components
          vtkFloatArray *newfa = vtkFloatArray::New();
          newfa->DeepCopy(fa);
          this->Superclass::VtkObjectPtr = newfa;
          break;
          }
        }
#endif
    case LABELLIST:
    case SCALARLIST:
    case STRINGLIST:
      this->Superclass::VtkObjectPtr = value.ToVTKObject();
      this->Superclass::VtkObjectPtr->Register(0);
      break;
    case LABELLISTLIST:
      this->LabelListListPtr = new vtkNewFoamIntVectorVector(*value.LabelListListPtr);
      break;
    case ENTRYVALUELIST:
      {
      const size_t nValues = value.EntryValuePtrs->size();
      this->EntryValuePtrs = new vtkstd::vector<vtkNewFoamEntryValue*>(nValues);
      for (size_t valueI = 0; valueI < nValues; valueI++)
        {
        this->EntryValuePtrs->operator[](valueI) = new vtkNewFoamEntryValue(
            *value.EntryValuePtrs->operator[](valueI), upperEntryPtr);
        }
      }
      break;
    case DICTIONARY:
      // UpperEntryPtr is null when called from vtkNewFoamDict constructor
      if (this->UpperEntryPtr != NULL)
        {
        this->DictPtr = new vtkNewFoamDict(*value.DictPtr,
            this->UpperEntryPtr->GetUpperDictPtr());
        }
      else
        {
        this->DictPtr = NULL;
        }
      break;
    case EMPTYLIST:
      break;
      // required to suppress the 'enumeration value not handled' warning by
      // g++ when compiled with -Wall
    default:
      break;
    }
}

void vtkNewFoamEntryValue::Clear()
{
  this->IsUniform = UNDEFINED;
  if (this->Managed)
    {
    switch (this->Superclass::Type)
      {
      case LABELLIST:
      case SCALARLIST:
      case VECTORLIST:
      case STRINGLIST:
        this->VtkObjectPtr->Delete();
        break;
      case LABELLISTLIST:
        delete this->LabelListListPtr;
        break;
      case ENTRYVALUELIST:
        for (size_t valueI = 0; valueI < this->EntryValuePtrs->size() ; valueI++)
          {
          delete this->EntryValuePtrs->operator[](valueI);
          }
        delete this->EntryValuePtrs;
        break;
      case DICTIONARY:
        delete this->DictPtr;
        break;
        // required to suppress the 'enumeration value not handled' warning by
        // g++ when compiled with -Wall
      default:
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
void vtkNewFoamEntryValue::ReadList(vtkNewFoamIOobject& io)
{
  vtkNewFoamToken currToken;
  io.Read(currToken);

  // initial guess of the list type
  if (currToken.GetType() == this->Superclass::LABEL)
    {
    // if the first token is of type LABEL it might be either an element of
    // a labelList or the size of a sublist so proceed to the next token
    vtkNewFoamToken nextToken;
    if (!io.Read(nextToken))
      {
      throw vtkNewFoamError() << "Unexpected EOF";
      }
    if (nextToken.GetType() == this->Superclass::LABEL)
      {
      this->Superclass::LabelListPtr = vtkIntArray::New();
      this->Superclass::LabelListPtr->InsertNextValue(currToken.To<int>());
      this->Superclass::LabelListPtr->InsertNextValue(nextToken.To<int>());
      this->Superclass::Type = LABELLIST;
      }
    else if (nextToken.GetType() == this->Superclass::SCALAR)
      {
      this->Superclass::ScalarListPtr = vtkFloatArray::New();
      this->Superclass::ScalarListPtr->InsertNextValue(currToken.To<float>());
      this->Superclass::ScalarListPtr->InsertNextValue(nextToken.To<float>());
      this->Superclass::Type = SCALARLIST;
      }
    else if (nextToken == '(') // list of list: read recursively
      {
      this->Superclass::EntryValuePtrs = new vtkstd::vector<vtkNewFoamEntryValue*>;
      this->Superclass::EntryValuePtrs->push_back(new vtkNewFoamEntryValue(
          this->UpperEntryPtr));
      this->Superclass::EntryValuePtrs->back()->ReadList(io);
      this->Superclass::Type = ENTRYVALUELIST;
      }
    else if (nextToken == ')') // list with only one label element
      {
      this->Superclass::LabelListPtr = vtkIntArray::New();
      this->Superclass::LabelListPtr->SetNumberOfValues(1);
      this->Superclass::LabelListPtr->SetValue(0, currToken.To<int>());
      this->Superclass::Type = LABELLIST;
      return;
      }
    else
      {
      throw vtkNewFoamError() << "Expected number, '(' or ')', found "
      << nextToken;
      }
    }
  else if (currToken.GetType() == this->Superclass::SCALAR)
    {
    this->Superclass::ScalarListPtr = vtkFloatArray::New();
    this->Superclass::ScalarListPtr->InsertNextValue(currToken.To<float>());
    this->Superclass::Type = SCALARLIST;
    }
  // if the first word is a string we have to read another token to determine
  // if the first word is a keyword for the following dictionary
  else if (currToken.IsWordOrString())
    {
    vtkNewFoamToken nextToken;
    if (!io.Read(nextToken))
      {
      throw vtkNewFoamError() << "Unexpected EOF";
      }
    if (nextToken.IsWordOrString()) // list of strings
      {
      this->Superclass::StringListPtr = vtkStringArray::New();
      this->Superclass::StringListPtr->InsertNextValue(currToken.ToStdString());
      this->Superclass::StringListPtr->InsertNextValue(nextToken.ToStdString());
      this->Superclass::Type = STRINGLIST;
      }
    // dictionary with the already read word or string as the first keyword
    else if (nextToken == '{')
      {
      this->ReadDictionary(io, currToken);
      // the dictionary read as list has the entry terminator ';' so
      // we have to skip it
      return;
      }
    else if (nextToken == ')') // list with only one string element
      {
      this->Superclass::StringListPtr = vtkStringArray::New();
      this->Superclass::StringListPtr->SetNumberOfValues(1);
      this->Superclass::StringListPtr->SetValue(0, currToken.ToStdString());
      this->Superclass::Type = STRINGLIST;
      return;
      }
    else
      {
      throw vtkNewFoamError() << "Expected string, '{' or ')', found "
      << nextToken;
      }
    }
  // list of lists or dictionaries: read recursively
  else if (currToken == '(' || currToken == '{')
    {
    this->Superclass::EntryValuePtrs = new vtkstd::vector<vtkNewFoamEntryValue*>;
    this->Superclass::EntryValuePtrs->push_back(new vtkNewFoamEntryValue(
        this->UpperEntryPtr));
    if(currToken == '(')
      {
      this->Superclass::EntryValuePtrs->back()->ReadList(io);
      }
    else // currToken == '{'
      {
      this->Superclass::EntryValuePtrs->back()->ReadDictionary(io, vtkNewFoamToken());
      }
    // read all the following values as arbitrary entryValues
    // the alphaContactAngle b.c. in multiphaseInterFoam/damBreak4phase
    // reaquires this treatment (reading by readList() is not enough)
    do
      {
      this->Superclass::EntryValuePtrs->push_back(new vtkNewFoamEntryValue(
          this->UpperEntryPtr));
      this->Superclass::EntryValuePtrs->back()->Read(io);
      } while (*this->Superclass::EntryValuePtrs->back() != ')'
        && *this->Superclass::EntryValuePtrs->back() != '}'
        && *this->Superclass::EntryValuePtrs->back() != ';');

    if (*this->Superclass::EntryValuePtrs->back() != ')')
      {
      throw vtkNewFoamError() << "Expected ')' before "
          << *this->Superclass::EntryValuePtrs->back();
      }

    // delete ')'
    delete this->Superclass::EntryValuePtrs->back();
    this->EntryValuePtrs->pop_back();
    this->Superclass::Type = ENTRYVALUELIST;
    return;
    }
  else if (currToken == ')') // empty list
    {
    this->Superclass::Type = EMPTYLIST;
    return;
    }
  // FIXME: may (or may not) need identifier handling

  while (io.Read(currToken) && currToken != ')')
    {
    if (this->Superclass::Type == LABELLIST)
      {
      if (currToken.GetType() == this->Superclass::SCALAR)
        {
        // switch to scalarList
        // LabelListPtr and ScalarListPtr are packed into a single union so
        // we need a temprary pointer
        vtkFloatArray* slPtr = vtkFloatArray::New();
        const int size = this->Superclass::LabelListPtr->GetNumberOfTuples();
        slPtr->SetNumberOfValues(size + 1);
        for (int i = 0; i < size; i++)
          {
          slPtr->SetValue(i,
              static_cast<float>(this->Superclass::LabelListPtr->GetValue(i)));
          }
        this->LabelListPtr->Delete();
        slPtr->SetValue(size, currToken.To<float>());
        // copy after LabelListPtr is deleted
        this->Superclass::ScalarListPtr = slPtr;
        this->Superclass::Type = SCALARLIST;
        }
      else if (currToken.GetType() == this->Superclass::LABEL)
        {
        this->Superclass::LabelListPtr->InsertNextValue(currToken.To<int>());
        }
      else
        {
        throw vtkNewFoamError() << "Expected a number, found " << currToken;
        }
      }
    else if (this->Superclass::Type == this->Superclass::SCALARLIST)
      {
      if (currToken.Is<float>())
        {
        this->Superclass::ScalarListPtr->InsertNextValue(currToken.To<float>());
        }
      else
        {
        throw vtkNewFoamError() << "Expected a number, found " << currToken;
        }
      }
    else if (this->Superclass::Type == this->Superclass::STRINGLIST)
      {
      if (currToken.IsWordOrString())
        {
        this->Superclass::StringListPtr->InsertNextValue(currToken.ToStdString());
        }
      else
        {
        throw vtkNewFoamError() << "Expected a string, found " << currToken;
        }
      }
    else if (this->Superclass::Type == this->Superclass::ENTRYVALUELIST)
      {
      if (currToken.GetType() == this->Superclass::LABEL)
        {
        // skip the number of elements to make things simple
        if (!io.Read(currToken))
          {
          throw vtkNewFoamError() << "Unexpected EOF";
          }
        }
      if (currToken != '(')
        {
        throw vtkNewFoamError() << "Expected '(', found " << currToken;
        }
      this->Superclass::EntryValuePtrs->push_back(new vtkNewFoamEntryValue(this->UpperEntryPtr));
      this->Superclass::EntryValuePtrs->back()->ReadList(io);
      }
    else
      {
      throw vtkNewFoamError() << "Unexpected token " << currToken;
      }
    }

  if (this->Superclass::Type == this->Superclass::LABELLIST)
    {
    this->Superclass::LabelListPtr->Squeeze();
    }
  else if (this->Superclass::Type == this->Superclass::SCALARLIST)
    {
    this->Superclass::ScalarListPtr->Squeeze();
    }
  else if (this->Superclass::Type == this->Superclass::STRINGLIST)
    {
    this->Superclass::StringListPtr->Squeeze();
    }
}

// a list of dictionaries is actually read as a dictionary
void vtkNewFoamEntryValue::ReadDictionary(vtkNewFoamIOobject& io,
    const vtkNewFoamToken& firstKeyword)
{
  this->Superclass::DictPtr = new vtkNewFoamDict(this->UpperEntryPtr->GetUpperDictPtr());
  this->Superclass::Type = this->Superclass::DICTIONARY;
  this->Superclass::DictPtr->Read(io, true, firstKeyword);
}

// guess the type of the given entry value and read it
void vtkNewFoamEntryValue::Read(vtkNewFoamIOobject& io)
{
  vtkNewFoamToken currToken;
  if (!io.Read(currToken))
    {
    throw vtkNewFoamError() << "Unexpected EOF";
    }

  if (currToken == '{')
    {
    this->ReadDictionary(io, vtkNewFoamToken());
    return;
    }
  // for reading sublist from vtkNewFoamEntryValue::readList() or there
  // are cases where lists without the (non)uniform keyword appear
  // (e. g. coodles/pitsDaily/0/U, uniformFixedValue b.c.)
  else if (currToken == '(')
    {
    this->ReadList(io);
    return;
    }
  else if (currToken == '[')
    {
    this->ReadDimensionSet(io);
    return;
    }
  else if (currToken == "uniform")
    {
    if (!io.Read(currToken))
      {
      throw vtkNewFoamError()
      << "Expected a uniform value or a list, found unexpected EOF";
      }
    if (currToken == '(')
      {
      this->ReadList(io);
      }
    else if (currToken.GetType() == this->Superclass::LABEL
        || currToken.GetType() == this->Superclass::SCALAR
        || currToken.IsWordOrString()
        || currToken.GetType() == this->Superclass::IDENTIFIER)
      {
      this->Superclass::operator=(currToken);
      }
    else // unexpected punctuation token
      {
      throw vtkNewFoamError() << "Expected number, string or (, found "
      << currToken;
      }
    this->IsUniform = UNIFORM;
    }
  else if (currToken == "nonuniform")
    {
    if (!io.Read(currToken))
      {
      throw vtkNewFoamError() << "Expected list type specifier, found EOF";
      }
    this->IsUniform = NONUNIFORM;
    if (currToken == "List<scalar>")
      {
      this->ReadNonuniformList<SCALARLIST, listTraits<vtkFloatArray, float> >(io);
      }
    else if (currToken == "List<sphericalTensor>")
      {
      this->ReadNonuniformList<VECTORLIST,
      vectorListTraits<vtkFloatArray, float, 1, false> >(io);
      }
    else if (currToken == "List<vector>")
      {
      this->ReadNonuniformList<VECTORLIST,
      vectorListTraits<vtkFloatArray, float, 3, false> >(io);
      }
    else if (currToken == "List<symmTensor>")
      {
      this->ReadNonuniformList<VECTORLIST,
      vectorListTraits<vtkFloatArray, float, 6, false> >(io);
      }
    else if (currToken == "List<tensor>")
      {
      this->ReadNonuniformList<VECTORLIST,
      vectorListTraits<vtkFloatArray, float, 9, false> >(io);
      }
    // List<bool> is read as List<label>
    else if (currToken == "List<label>" || currToken == "List<bool>")
      {
      this->ReadNonuniformList<LABELLIST, listTraits<vtkIntArray, int> >(io);
      }
    // an empty list doesn't have a list type specifier
    else if (currToken.GetType() == this->Superclass::LABEL
        && currToken.To<int>() == 0)
      {
      this->Superclass::Type = this->Superclass::EMPTYLIST;
      if(io.GetFormat() == vtkNewFoamIOobject::ASCII)
        {
        io.ReadExpecting('(');
        io.ReadExpecting(')');
        }
      }
    else if (currToken.GetType() == this->Superclass::IDENTIFIER)
      {
      this->Superclass::operator=(currToken);
      }
    else
      {
      throw vtkNewFoamError() << "Unsupported nonuniform list type " << currToken;
      }
    }
  // the followings read nonuniform lists without setting isUniform to
  // NONUNIFORM
  else if (currToken == "List<scalar>")
    {
    this->ReadNonuniformList<SCALARLIST, listTraits<vtkFloatArray, float> >(io);
    }
  else if (currToken == "List<sphericalTensor>")
    {
    this->ReadNonuniformList<VECTORLIST,
      vectorListTraits<vtkFloatArray, float, 1, false> >(io);
    }
  else if (currToken == "List<vector>")
    {
    this->ReadNonuniformList<VECTORLIST,
      vectorListTraits<vtkFloatArray, float, 3, false> >(io);
    }
  else if (currToken == "List<symmTensor>")
    {
    this->ReadNonuniformList<VECTORLIST,
      vectorListTraits<vtkFloatArray, float, 6, false> >(io);
    }
  else if (currToken == "List<tensor>")
    {
    this->ReadNonuniformList<VECTORLIST,
      vectorListTraits<vtkFloatArray, float, 9, false> >(io);
    }
  // zones have list without a uniform/nonuniform keyword
  // List<bool> is read as List<label>
  // (e. g. flipMap entry in faceZones)
  else if (currToken == "List<label>" || currToken == "List<bool>")
    {
    this->ReadNonuniformList<LABELLIST, listTraits<vtkIntArray, int> >(io);
    }
  else if (currToken.GetType() == this->Superclass::PUNCTUATION
      || currToken.GetType() == this->Superclass::LABEL || currToken.GetType()
      == this->Superclass::SCALAR || currToken.IsWordOrString()
      || currToken.GetType() == this->Superclass::IDENTIFIER)
    {
    this->Superclass::operator=(currToken);
    }
}

// read values of an entry
void vtkNewFoamEntry::Read(vtkNewFoamIOobject& io)
{
  for (;;)
    {
    this->Superclass::push_back(new vtkNewFoamEntryValue(this));
    this->Superclass::back()->Read(io);

    if (this->Superclass::size() >= 2)
      {
      vtkNewFoamEntryValue& secondLastValue =
          *this->Superclass::operator[](this->Superclass::size() - 2);
      if (secondLastValue.GetType() == vtkNewFoamToken::LABEL)
        {
        vtkNewFoamEntryValue& lastValue = *this->Superclass::back();

        // a zero-sized nonuniform list without prefixing "nonuniform"
        // keyword nor list type specifier (i. e. `0()';
        // e. g. simpleEngine/0/polyMesh/pointZones) requires special
        // care (one with nonuniform prefix is treated within
        // vtkNewFoamEntryValue::read()). still this causes errornous
        // behavior for `0 nonuniform 0()' but this should be extremely
        // rare
        if (lastValue.GetType() == vtkNewFoamToken::EMPTYLIST && secondLastValue
            == 0)
          {
          delete this->Superclass::back();
          this->Superclass::pop_back(); // delete the last value
          // mark new last value as empty
          this->Superclass::back()->SetEmptyList();
          }
        // for an exceptional expression of `LABEL{LABELorSCALAR}' without
        // type prefix (e. g. `2{-0}' in mixedRhoE B.C. in
        // rhopSonicFoam/shockTube)
        else if (lastValue.GetType() == vtkNewFoamToken::DICTIONARY)
          {
          if (lastValue.Dictionary().GetType() == vtkNewFoamToken::LABEL)
            {
            const int asize = secondLastValue.To<int>();
            // not using templated To<int>() for workarounding an error
            // on SunOS-CC5.6-dbg
            const int value = lastValue.Dictionary().GetToken().ToInt();
            // delete last two values
            delete this->Superclass::back();
            this->Superclass::pop_back();
            delete this->Superclass::back();
            this->Superclass::pop_back();
            // make new labelList
            this->Superclass::push_back(new vtkNewFoamEntryValue(this));
            this->Superclass::back()->MakeLabelList(value, asize);
            }
          else if (lastValue.Dictionary().GetType() == vtkNewFoamToken::SCALAR)
            {
            const int asize = secondLastValue.To<int>();
            // not using templated To<float>() for workarounding an error
            // on SunOS-CC5.6-dbg
            const float value = lastValue.Dictionary().GetToken().ToFloat();
            // delete last two values
            delete this->Superclass::back();
            this->Superclass::pop_back();
            delete this->Superclass::back();
            this->Superclass::pop_back();
            // make new labelList
            this->Superclass::push_back(new vtkNewFoamEntryValue(this));
            this->Superclass::back()->MakeScalarList(value, asize);
            }
          }
        }
      }

    if (this->Superclass::back()->GetType() == vtkNewFoamToken::IDENTIFIER)
      {
      // substitute identifier
      const vtkStdString identifier(this->Superclass::back()->ToStdString());
      vtkNewFoamEntryValue::uniformTypes isUniform
          = this->Superclass::back()->GetIsUniform();
      delete this->Superclass::back();
      this->Superclass::pop_back();

      for (const vtkNewFoamDict *uDictPtr = this->UpperDictPtr; uDictPtr != NULL;)
        {
        const vtkNewFoamEntry *identifiedEntry = uDictPtr->LookupLiteral(identifier);

        if (identifiedEntry != NULL)
          {
          if (isUniform != vtkNewFoamEntryValue::UNDEFINED
              && identifiedEntry->size() > 0
              && identifiedEntry->operator[](0)->GetIsUniform()
              != vtkNewFoamEntryValue::UNDEFINED)
            {
            throw vtkNewFoamError() << "duplicated list type specifiers (uniform/"
              "nonuniform) in the first substituted entry value of identifier $"
              << identifier;
            }

          for (size_t valueI = 0; valueI < identifiedEntry->size(); valueI++)
            {
            if (identifiedEntry->operator[](valueI)->GetType()
                == vtkNewFoamToken::DICTIONARY)
              {
              // in fact there's no problem with substituting
              // dictionary in this parser. Only for compatibility
              // with the genuine parser
              throw vtkNewFoamError() << "dictionary cannot be substituted as an "
                  "entry value. Try {$" << identifier << ";} instead of $"
                  << identifier << ".";
              }
            this->Superclass::push_back(new vtkNewFoamEntryValue(
                *identifiedEntry->operator[](valueI), this));
            }

          if (isUniform != vtkNewFoamEntryValue::UNDEFINED)
            {
            if (identifiedEntry->size() > 0)
              {
              this->Superclass::operator[](this->Superclass::size()
                  - identifiedEntry->size())->SetIsUniform(isUniform);
              }
            else
              {
              this->Superclass::push_back(new vtkNewFoamEntryValue(this));
              this->Superclass::back()->SetIsUniform(isUniform);
              }
            }
          break;
          }
        else
          {
          uDictPtr = uDictPtr->GetUpperDictPtr();
          if (uDictPtr == NULL)
            {
            // if the substituting entry is not found append the identifier as
            // word
            this->Superclass::push_back(new vtkNewFoamEntryValue(this));
            this->Superclass::back()->SetWord("$" + identifier);
            this->Superclass::back()->SetIsUniform(isUniform);
            }
          }
        }
      }
    else if (*this->Superclass::back() == ';')
      {
      delete this->Superclass::back();
      this->Superclass::pop_back();
      break;
      }
    else if (this->Superclass::back()->GetType() == vtkNewFoamToken::DICTIONARY)
      {
      // subdictionary is not suffixed by an entry terminator ';'
      break;
      }
    else if (*this->Superclass::back() == '}' || *this->Superclass::back()
        == ')')
      {
      throw vtkNewFoamError() << "Unmatched " << *this->Superclass::back();
      }
    }
}

//-----------------------------------------------------------------------------
// vtkNewOpenFOAMReaderPrivate constructor and destructor
vtkNewOpenFOAMReaderPrivate::vtkNewOpenFOAMReaderPrivate()
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
#if 0
  this->ReciprocalDelta = NULL;
#endif
  this->FaceOwner = NULL;
  this->PointZoneMesh = NULL;
  this->FaceZoneMesh = NULL;
  this->CellZoneMesh = NULL;

  // for decomposing polyhedra
  this->NumAdditionalCells = 0;
  this->AdditionalCellIds = NULL;
  this->NumAdditionalCells = NULL;
  this->AdditionalCellPoints = NULL;
}

vtkNewOpenFOAMReaderPrivate::~vtkNewOpenFOAMReaderPrivate()
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

void vtkNewOpenFOAMReaderPrivate::PrintSelf(ostream &os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
  os << indent << "Case Path: "
      << (this->CasePath.length() ? this->CasePath.c_str() : "(none)") << endl;
  os << indent << "Region Name: "
      << (this->RegionName.length() ? this->RegionName.c_str() : "(none)")
      << endl;
  os << indent << "Processor Name: "
      << (this->ProcessorName.length() ? this->ProcessorName.c_str() : "(none)")
      << endl;
  if (this->TimeValues)
    {
    os << indent << "Number of Time Steps: "
        << this->TimeValues->GetNumberOfTuples() << endl;
    }
  os << indent << "Time Step: " << this->TimeStep << endl;
  os << indent << "Number of Cells: " << this->NumCells << endl;
  os << indent << "Number of Points: " << this->NumPoints << endl;
}

void vtkNewOpenFOAMReaderPrivate::ClearInternalMeshes()
{
  if (this->FaceOwner != NULL)
    {
    this->FaceOwner->Delete();
    this->FaceOwner = NULL;
    }
  if (this->InternalMesh != NULL)
    {
    this->InternalMesh->Delete();
    this->InternalMesh = NULL;
    }
  if (this->AdditionalCellIds != NULL)
    {
    this->AdditionalCellIds->Delete();
    this->AdditionalCellIds = NULL;
    }
  if (this->NumAdditionalCells != NULL)
    {
    this->NumAdditionalCells->Delete();
    this->NumAdditionalCells = NULL;
    }
  delete this->AdditionalCellPoints;
  this->AdditionalCellPoints = NULL;

  if (this->PointZoneMesh != NULL)
    {
    this->PointZoneMesh->Delete();
    this->PointZoneMesh = NULL;
    }
  if (this->FaceZoneMesh != NULL)
    {
    this->FaceZoneMesh->Delete();
    this->FaceZoneMesh = NULL;
    }
  if (this->CellZoneMesh != NULL)
    {
    this->CellZoneMesh->Delete();
    this->CellZoneMesh = NULL;
    }
}

void vtkNewOpenFOAMReaderPrivate::ClearBoundaryMeshes()
{
  if (this->BoundaryMesh != NULL)
    {
    this->BoundaryMesh->Delete();
    this->BoundaryMesh = NULL;
    }

  delete this->BoundaryPointMap;
  this->BoundaryPointMap = NULL;
#if 0
  delete this->ReciprocalDelta;
  this->ReciprocalDelta = NULL;
#endif

  if (this->InternalPoints != NULL)
    {
    this->InternalPoints->Delete();
    this->InternalPoints = NULL;
    }
  if (this->AllBoundaries != NULL)
    {
    this->AllBoundaries->Delete();
    this->AllBoundaries = NULL;
    }
  if (this->AllBoundariesPointMap != NULL)
    {
    this->AllBoundariesPointMap->Delete();
    this->AllBoundariesPointMap = NULL;
    }
}

void vtkNewOpenFOAMReaderPrivate::ClearMeshes()
{
  this->ClearInternalMeshes();
  this->ClearBoundaryMeshes();
}

void vtkNewOpenFOAMReaderPrivate::SetTimeValue(const double requestedTime)
{
  const int nTimeValues = this->TimeValues->GetNumberOfTuples();
  if (nTimeValues > 0)
    {
    int minTimeI = 0;
    double minTimeDiff = fabs(this->TimeValues->GetValue(0) - requestedTime);
    for (int timeI = 1; timeI < nTimeValues; timeI++)
      {
      const double timeDiff(fabs(this->TimeValues->GetValue(timeI)
          - requestedTime));
      if (timeDiff < minTimeDiff)
        {
        minTimeI = timeI;
        minTimeDiff = timeDiff;
        }
      }
    this->SetTimeStep(minTimeI); // set Modified() if TimeStep changed
    }
}

//-----------------------------------------------------------------------------
void vtkNewOpenFOAMReaderPrivate::SetupInformation(const vtkStdString &casePath,
    const vtkStdString &regionName, const vtkStdString &procName,
    vtkNewOpenFOAMReaderPrivate *master)
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
void vtkNewOpenFOAMReaderPrivate::GetFieldNames(const vtkStdString &tempPath,
    const bool isLagrangian, vtkStringArray *cellObjectNames,
    vtkStringArray *pointObjectNames)
{
  // open the directory and get num of files
  vtkDirectory *directory = vtkDirectory::New();
  if (!directory->Open(tempPath.c_str()))
    {
    // no data
    directory->Delete();
    return;
    }

  // loop over all files and locate valid fields
  int nFieldFiles = directory->GetNumberOfFiles();
  for (int j = 0; j < nFieldFiles; j++)
    {
    const vtkStdString fieldFile(directory->GetFile(j));
    const size_t len = fieldFile.length();

    // excluded extensions cf. src/OpenFOAM/OSspecific/Unix/Unix.C
    if (!directory->FileIsDirectory(fieldFile.c_str()) && fieldFile.substr(len
        - 1) != "~" && (len < 4 || (fieldFile.substr(len - 4) != ".bak"
        && fieldFile.substr(len - 4) != ".BAK" && fieldFile.substr(len - 4)
        != ".old")) && (len < 5 || fieldFile.substr(len - 5) != ".save"))
      {
      vtkNewFoamIOobject io(this->CasePath,
          this->Parent->GetIsSinglePrecisionBinary() != 0);
      if (io.Open(tempPath + "/" + fieldFile)) // file exists and readable
        {
        const vtkStdString& cn = io.GetClassName();
        if (isLagrangian)
          {
          if (cn == "labelField" || cn == "scalarField" || cn == "vectorField"
              || cn == "sphericalTensorField" || cn == "symmTensorField" || cn
              == "tensorField")
            {
            // real file name
            this->LagrangianFieldFiles->InsertNextValue(fieldFile);
            // object name
            pointObjectNames->InsertNextValue(io.GetObjectName());
            }
          }
        else
          {
          if (cn == "volScalarField" || cn == "pointScalarField" || cn
              == "volVectorField" || cn == "pointVectorField" || cn
              == "volSphericalTensorField" || cn == "pointSphericalTensorField"
              || cn == "volSymmTensorField" || cn == "pointSymmTensorField"
              || cn == "volTensorField" || cn == "pointTensorField")
            {
            if (cn.substr(0, 3) == "vol")
              {
              // real file name
              this->VolFieldFiles->InsertNextValue(fieldFile);
              // object name
              cellObjectNames->InsertNextValue(io.GetObjectName());
              }
            else
              {
              this->PointFieldFiles->InsertNextValue(fieldFile);
              pointObjectNames->InsertNextValue(io.GetObjectName());
              }
            }
          }
        io.Close();
        }
      }
    }
  // inserted objects are squeezed later in SortFieldFiles()
  directory->Delete();
}

//-----------------------------------------------------------------------------
// locate laglangian clouds
void vtkNewOpenFOAMReaderPrivate::LocateLagrangianClouds(
    vtkStringArray *lagrangianObjectNames, const vtkStdString &timePath)
{
  vtkDirectory *directory = vtkDirectory::New();
  if (directory->Open((timePath + this->RegionPath() + "/lagrangian").c_str()))
    {
    // search for sub-clouds (OF 1.5 format)
    const int nFiles = directory->GetNumberOfFiles();
    bool isSubCloud = false;
    for (int fileI = 0; fileI < nFiles; fileI++)
      {
      const vtkStdString fileNameI(directory->GetFile(fileI));
      if (fileNameI != "." && fileNameI != ".."
          && directory->FileIsDirectory(fileNameI.c_str()))
        {
        vtkNewFoamIOobject io(this->CasePath,
            this->Parent->GetIsSinglePrecisionBinary() != 0);
        const vtkStdString subCloudName(this->RegionPrefix() + "lagrangian/"
            + fileNameI);
        const vtkStdString subCloudFullPath(timePath + "/" + subCloudName);
        // lagrangian positions. there are many concrete class names
        // e. g. Cloud<parcel>, basicKinematicCloud etc.
        if ((io.Open(subCloudFullPath + "/positions")
            || io.Open(subCloudFullPath + "/positions.gz")) && io.GetClassName().find("Cloud") != vtkStdString::npos && io.GetObjectName()
            == "positions")
          {
          isSubCloud = true;
          // a lagrangianPath has to be in a bit different format from
          // subCloudName to make the "lagrangian" reserved path
          // component and a mesh region with the same name
          // distinguishable later
          const vtkStdString subCloudPath(this->RegionName + "/lagrangian/"
              + fileNameI);
          if (this->Parent->LagrangianPaths->LookupValue(subCloudPath) == -1)
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
    if (!isSubCloud)
      {
      vtkNewFoamIOobject io(this->CasePath,
          this->Parent->GetIsSinglePrecisionBinary() != 0);
      const vtkStdString cloudName(this->RegionPrefix() + "lagrangian");
      const vtkStdString cloudFullPath(timePath + "/" + cloudName);
      if ((io.Open(cloudFullPath + "/positions") || io.Open(cloudFullPath
          + "/positions.gz")) && io.GetClassName().find("Cloud") != vtkStdString::npos && io.GetObjectName()
          == "positions")
        {
        const vtkStdString cloudPath(this->RegionName + "/lagrangian");
        if (this->Parent->LagrangianPaths->LookupValue(cloudPath) == -1)
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
void vtkNewOpenFOAMReaderPrivate::SortFieldFiles(vtkStringArray *selections,
    vtkStringArray *files, vtkStringArray *objects)
{
  objects->Squeeze();
  files->Squeeze();
  vtkSortDataArray::Sort(objects, files);
  for (int nameI = 0; nameI < objects->GetNumberOfValues(); nameI++)
    {
    selections->InsertNextValue(objects->GetValue(nameI));
    }
  objects->Delete();
}

//-----------------------------------------------------------------------------
// create field data lists and cell/point array selection lists
int vtkNewOpenFOAMReaderPrivate::MakeMetaDataAtTimeStep(
    vtkStringArray *cellSelectionNames, vtkStringArray *pointSelectionNames,
    vtkStringArray *lagrangianSelectionNames, const bool listNextTimeStep)
{
  // Read the patches from the boundary file into selection array
  if (this->PolyMeshFacesDir->GetValue(this->TimeStep)
      != this->BoundaryDict.TimeDir
      || this->Parent->PatchDataArraySelection->GetMTime()
          != this->Parent->PatchSelectionMTimeOld)
    {
    this->BoundaryDict.clear();
    this->BoundaryDict.TimeDir
        = this->PolyMeshFacesDir->GetValue(this->TimeStep);

    vtkNewFoamDict *boundaryDict = this->GatherBlocks("boundary", this->TimeStep);
    if (boundaryDict == NULL && listNextTimeStep
        && this->TimeValues->GetNumberOfTuples() >= 2 && this->TimeStep == 0)
      {
      this->BoundaryDict.TimeDir
          = this->PolyMeshFacesDir->GetValue(1);
      boundaryDict = this->GatherBlocks("boundary", 1);
      }
    if (boundaryDict != NULL)
      {
      // Add the internal mesh by default always
      const vtkStdString
          internalMeshName(this->RegionPrefix() + "internalMesh");
      if (this->RegionName == ""
          && this->Parent->Readers->GetNumberOfItems() > 1
          && !this->Parent->PatchDataArraySelection
          ->ArrayExists(internalMeshName.c_str()))
        {
        this->Parent->PatchDataArraySelection->DisableArray(
            internalMeshName.c_str());
        }
      else
        {
        this->Parent->PatchDataArraySelection
            ->AddArray(internalMeshName.c_str());
        }
      this->InternalMeshSelectionStatus
          = this->Parent->GetPatchArrayStatus(internalMeshName.c_str());

      // iterate through each entry in the boundary file
      int allBoundariesNextStartFace = 0;
      this->BoundaryDict.resize(boundaryDict->size());
      for (size_t i = 0; i < boundaryDict->size(); i++)
        {
        vtkNewFoamEntry *boundaryEntryI = boundaryDict->operator[](i);
        const vtkNewFoamEntry *nFacesEntry = boundaryEntryI->Dictionary().Lookup("nFaces");
        if (nFacesEntry == NULL)
          {
          vtkErrorMacro(<< "nFaces entry not found in boundary entry "
              << boundaryEntryI->GetKeyword().c_str());
          delete boundaryDict;
          return 0;
          }
        const int nFaces = nFacesEntry->ToInt();

        // extract name of the current patch for insertion
        const vtkStdString &boundaryNameI = boundaryEntryI->GetKeyword();

        // create BoundaryDict entry
        vtkNewFoamBoundaryEntry &BoundaryEntryI = this->BoundaryDict[i];
        BoundaryEntryI.NFaces = nFaces;
        BoundaryEntryI.BoundaryName = boundaryNameI;
        const vtkNewFoamEntry *startFaceEntry = boundaryEntryI->Dictionary().Lookup("startFace");
        if (startFaceEntry == NULL)
          {
          vtkErrorMacro(<< "startFace entry not found in boundary entry "
              << boundaryEntryI->GetKeyword().c_str());
          delete boundaryDict;
          return 0;
          }
        BoundaryEntryI.StartFace = startFaceEntry->ToInt();
        const vtkNewFoamEntry *typeEntry = boundaryEntryI->Dictionary().Lookup("type");
        if (typeEntry == NULL)
          {
          vtkErrorMacro(<< "type entry not found in boundary entry "
              << boundaryEntryI->GetKeyword().c_str());
          delete boundaryDict;
          return 0;
          }
        BoundaryEntryI.AllBoundariesStartFace = allBoundariesNextStartFace;
        const vtkStdString typeNameI(typeEntry->ToWord());
        // if the basic type of the patch is one of the followings the
        // point-filtered values at patches are overridden by patch values
        if (typeNameI == "patch" || typeNameI == "wall")
          {
          BoundaryEntryI.BoundaryType = vtkNewFoamBoundaryEntry::PHYSICAL;
          allBoundariesNextStartFace += nFaces;
          }
        else if (typeNameI == "processor")
          {
          BoundaryEntryI.BoundaryType = vtkNewFoamBoundaryEntry::PROCESSOR;
          allBoundariesNextStartFace += nFaces;
          }
        else
          {
          BoundaryEntryI.BoundaryType = vtkNewFoamBoundaryEntry::GEOMETRICAL;
          }
        BoundaryEntryI.IsActive = false;

        // always hide processor patches for decomposed cases to keep
        // vtkAppendCompositeDataLeaves happy
        if (this->ProcessorName != "" && BoundaryEntryI.BoundaryType
            == vtkNewFoamBoundaryEntry::PROCESSOR)
          {
          continue;
          }
        const vtkStdString selectionName(this->RegionPrefix() + boundaryNameI);
        if (this->Parent->PatchDataArraySelection->
        ArrayExists(selectionName.c_str()))
          {
          // Mark boundary if selected for display
          if (this->Parent->GetPatchArrayStatus(selectionName.c_str()))
            {
            BoundaryEntryI.IsActive = true;
            }
          }
        else
          {
          // add patch to list with selection status turned off:
          // the patch is added to list even if its size is zero
          this->Parent->PatchDataArraySelection->DisableArray(selectionName.c_str());
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
  if (listNextTimeStep)
    {
    this->Parent->LagrangianPaths->Initialize();
    }
  vtkStringArray *lagrangianObjectNames = vtkStringArray::New();
  this->LocateLagrangianClouds(lagrangianObjectNames, timePath);

  // if the requested timestep is 0 then we also look at the next
  // timestep to add extra objects that don't exist at timestep 0 into
  // selection lists. Note the ObjectNames array will be recreated in
  // RequestData() so we don't have to worry about duplicated fields.
  if (listNextTimeStep && this->TimeValues->GetNumberOfTuples() >= 2
      && this->TimeStep == 0)
    {
    const vtkStdString timePath2(this->TimePath(1));
    this->GetFieldNames(timePath2 + this->RegionPath(), false, cellObjectNames,
        pointObjectNames);
    // if lagrangian clouds were not found at timestep 0
    if (this->Parent->LagrangianPaths->GetNumberOfTuples() == 0)
      {
      this->LocateLagrangianClouds(lagrangianObjectNames, timePath2);
      }
    }

  // sort array names
  this->SortFieldFiles(cellSelectionNames, this->VolFieldFiles, cellObjectNames);
  this->SortFieldFiles(pointSelectionNames, this->PointFieldFiles,
      pointObjectNames);
  this->SortFieldFiles(lagrangianSelectionNames, this->LagrangianFieldFiles,
      lagrangianObjectNames);

  return 1;
}

//-----------------------------------------------------------------------------
// list time directories according to controlDict
bool vtkNewOpenFOAMReaderPrivate::ListTimeDirectoriesByControlDict(
    vtkNewFoamDict* dictPtr)
{
  vtkNewFoamDict& dict = *dictPtr;

  const vtkNewFoamEntry *startTimeEntry = dict.Lookup("startTime");
  if (startTimeEntry == NULL)
    {
    vtkErrorMacro(<< "startTime entry not found in controlDict");
    return false;
    }
  // using double to precisely handle time values
  const double startTime = startTimeEntry->ToDouble();

  const vtkNewFoamEntry *endTimeEntry = dict.Lookup("endTime");
  if (endTimeEntry == NULL)
    {
    vtkErrorMacro(<< "endTime entry not found in controlDict");
    return false;
    }
  const double endTime = endTimeEntry->ToDouble();

  const vtkNewFoamEntry *deltaTEntry = dict.Lookup("deltaT");
  if (deltaTEntry == NULL)
    {
    vtkErrorMacro(<< "deltaT entry not found in controlDict");
    return false;
    }
  const double deltaT = deltaTEntry->ToDouble();

  const vtkNewFoamEntry *writeIntervalEntry = dict.Lookup("writeInterval");
  if (writeIntervalEntry == NULL)
    {
    vtkErrorMacro(<< "writeInterval entry not found in controlDict");
    return false;
    }
  const double writeInterval = writeIntervalEntry->ToDouble();

  const vtkNewFoamEntry *timeFormatEntry = dict.Lookup("timeFormat");
  if (timeFormatEntry == NULL)
    {
    vtkErrorMacro(<< "timeFormat entry not found in controlDict");
    return false;
    }
  const vtkStdString timeFormat(timeFormatEntry->ToWord());

  const vtkNewFoamEntry *timePrecisionEntry = dict.Lookup("timePrecision");
  const int timePrecision // default is 6
      = (timePrecisionEntry != NULL ? timePrecisionEntry->ToInt() : 6);

  // calculate the time step increment based on type of run
  const vtkNewFoamEntry *writeControlEntry = dict.Lookup("writeControl");
  if (writeControlEntry == NULL)
    {
    vtkErrorMacro(<< "writeControl entry not found in controlDict");
    return false;
    }
  const vtkStdString writeControl(writeControlEntry->ToWord());
  double timeStepIncrement;
  if (writeControl == "timeStep")
    {
    timeStepIncrement = writeInterval * deltaT;
    }
  else if (writeControl == "runTime" || writeControl == "adjustableRunTime")
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
  if (timeFormat == "general")
    {
    parser.setf(vtksys_ios::ios_base::fmtflags(0), vtksys_ios::ios_base::floatfield);
    }
  else if (timeFormat == "fixed")
    {
    parser.setf(vtksys_ios::ios_base::fmtflags(vtksys_ios::ios_base::fixed),
        vtksys_ios::ios_base::floatfield);
#ifdef _MSC_VER
    correctExponent = false;
#endif
    }
  else if (timeFormat == "scientific")
    {
    parser.setf(vtksys_ios::ios_base::fmtflags(vtksys_ios::ios_base::scientific),
        vtksys_ios::ios_base::floatfield);
    }
  else
    {
    vtkWarningMacro("Warning: unsupported time format. Assuming general.");
    parser.setf(vtksys_ios::ios_base::fmtflags(0), vtksys_ios::ios_base::floatfield);
    }
  parser.precision(timePrecision);

  for (int i = 0; i < tempNumTimeSteps; i++)
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
    if (test->Open((this->CasePath + parser.str()).c_str()))
      {
      this->TimeValues->InsertNextValue(tempStep);
      this->TimeNames->InsertNextValue(parser.str());
      }
    // necessary for reading the case/0 directory whatever the timeFormat is
    // based on Foam::Time::operator++() cf. src/OpenFOAM/db/Time/Time.C
    else if ((fabs(tempStep) < 1.0e-14L) // 10*SMALL
        && test->Open((this->CasePath + vtkStdString("0")).c_str()))
      {
      this->TimeValues->InsertNextValue(tempStep);
      this->TimeNames->InsertNextValue(vtkStdString("0"));
      }
    }
  test->Delete();
  this->TimeValues->Squeeze();
  this->TimeNames->Squeeze();

  if (this->TimeValues->GetNumberOfTuples() == 0)
    {
    // set the number of timesteps to 1 if the constant subdirectory exists
    test = vtkDirectory::New();
    if (test->Open((this->CasePath + "constant").c_str()))
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
bool vtkNewOpenFOAMReaderPrivate::ListTimeDirectoriesByInstances()
{
  // open the case directory
  vtkDirectory* test = vtkDirectory::New();
  if (!test->Open(this->CasePath.c_str()))
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
  for (int i = 0; i < nFiles; i++)
    {
    const vtkStdString dir = test->GetFile(i);
    if (test->FileIsDirectory(dir.c_str()))
      {
      // check if the name is convertible to a number
      bool isTimeDir = true;
      for (size_t j = 0; j < dir.length(); j++)
        {
        const char c = dir[j];
        if (!isdigit(c) && c != '+' && c != '-' && c != '.' && c != 'e' && c
            != 'E')
          {
          isTimeDir = false;
          break;
          }
        }
      if (!isTimeDir)
        {
        continue;
        }

      // convert to a number
#if VTK_FOAMFILE_LOCALE_WORKAROUND
      vtksys_ios::istringstream timeStream(dir);
      double timeValue;
      timeStream >> timeValue;

      // check if the value really was converted to a number
      if (timeStream.fail() || !timeStream.eof())
        {
        continue;
        }
#else
      char *endptr;
      double timeValue = strtod(dir.c_str(), &endptr);
      // check if the value really was converted to a number
      if (timeValue == 0.0 && endptr == dir.c_str())
        {
        continue;
        }
#endif

      // add to the instance list
      this->TimeValues->InsertNextValue(timeValue);
      this->TimeNames->InsertNextValue(dir);
      }
    }
  test->Delete();

  this->TimeValues->Squeeze();
  this->TimeNames->Squeeze();

  if (this->TimeValues->GetNumberOfTuples() > 1)
    {
    // sort the detected time directories
    vtkSortDataArray::Sort(this->TimeValues, this->TimeNames);

    // if there are duplicated timeValues found, remove duplicates
    // (e.g. "0" and "0.000")
    for (int timeI = 1; timeI < this->TimeValues->GetNumberOfTuples();)
      {
      // compare by exact match
      if (this->TimeValues->GetValue(timeI - 1)
          == this->TimeValues->GetValue(timeI))
        {
        vtkWarningMacro(<<"Different time directories with the same time value "
            << this->TimeNames->GetValue(timeI - 1).c_str() << " and "
            << this->TimeNames->GetValue(timeI).c_str() << " found. "
            << this->TimeNames->GetValue(timeI).c_str() << " will be ignored.");
        this->TimeValues->RemoveTuple(timeI);
        // vtkStringArray does not have RemoveTuple()
        for (int timeJ = timeI + 1; timeJ
            < this->TimeNames->GetNumberOfTuples(); timeJ++)
          {
          this->TimeNames->SetValue(timeJ - 1, this->TimeNames->GetValue(timeJ));
          }
        this->TimeNames->Resize(this->TimeNames->GetNumberOfTuples() - 1);
        }
      else
        {
        timeI++;
        }
      }
    }
  else if (this->TimeValues->GetNumberOfTuples() == 0)
    {
    // set the number of timesteps to 1 if the constant subdirectory exists
    test = vtkDirectory::New();
    if (test->Open((this->CasePath + "constant").c_str()))
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
bool vtkNewOpenFOAMReaderPrivate::MakeInformationVector(
    const vtkStdString &casePath, const vtkStdString &controlDictPath,
    const vtkStdString &procName, vtkNewOpenFOAMReader *parent)
{
  this->CasePath = casePath;
  this->ProcessorName = procName;
  this->Parent = parent;

  // list timesteps (skip parsing controlDict entirely if
  // ListTimeStepsByControlDict is set to 0)
  bool ret = false; // tentatively set to false to suppress warning by older compilers
  if (this->Parent->GetListTimeStepsByControlDict())
    {
    vtkNewFoamIOobject io(this->CasePath,
        this->Parent->GetIsSinglePrecisionBinary() != 0);

    // open and check if controlDict is readable
    if (!io.Open(controlDictPath))
      {
      vtkErrorMacro(<<"Error opening " << io.GetFileName().c_str() << ": "
          << io.GetError().c_str());
      return false;
      }
    vtkNewFoamDict dict;
    if (!dict.Read(io))
      {
      vtkErrorMacro(<<"Error reading line " << io.GetLineNumber()
          << " of " << io.GetFileName().c_str() << ": " << io.GetError().c_str());
      return false;
      }
    if (dict.GetType() != vtkNewFoamToken::DICTIONARY)
      {
      vtkErrorMacro(<<"The file type of " << io.GetFileName().c_str()
          << " is not a dictionary");
      return false;
      }

    const vtkNewFoamEntry *writeControlEntry = dict.Lookup("writeControl");
    if (writeControlEntry == NULL)
      {
      vtkErrorMacro(<< "writeControl entry not found in "
          << io.GetFileName().c_str());
      return false;
      }
    const vtkStdString writeControl(writeControlEntry->ToWord());

    // empty if not found
    const vtkNewFoamEntry *adjustTimeStepEntry = dict.Lookup("adjustTimeStep");
    const vtkStdString
        adjustTimeStep = adjustTimeStepEntry == NULL ? vtkStdString()
            : adjustTimeStepEntry->ToWord();

    // list time directories according to controlDict if (adjustTimeStep
    // writeControl) == (off, timeStep) or (on, adjustableRunTime); list
    // by time instances in the case directory otherwise (different behavior
    // from paraFoam)
    // valid switching words cf. src/OpenFOAM/db/Switch/Switch.C
    if ((((adjustTimeStep == "off" || adjustTimeStep == "no" || adjustTimeStep
        == "n" || adjustTimeStep == "false" || adjustTimeStep == "")
        && writeControl == "timeStep") || ((adjustTimeStep == "on"
        || adjustTimeStep == "yes" || adjustTimeStep == "y" || adjustTimeStep
        == "true") && writeControl == "adjustableRunTime")))
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

  if (!ret)
    {
    return ret;
    }

  // does not seem to be required even if number of timesteps reduced
  // upon refresh since ParaView rewinds TimeStep to 0, but for precaution
  if (this->TimeValues->GetNumberOfTuples() > 0)
    {
    if (this->TimeStep >= this->TimeValues->GetNumberOfTuples())
      {
      this->SetTimeStep(this->TimeValues->GetNumberOfTuples() - 1);
      }
    }
  else
    {
    this->SetTimeStep(0);
    }

  this->PopulatePolyMeshDirArrays();
  return ret;
}

//-----------------------------------------------------------------------------
void vtkNewOpenFOAMReaderPrivate::AppendMeshDirToArray(
    vtkStringArray* polyMeshDir, const vtkStdString &path,
    const vtkStdString &fileName, const int timeI)
{
  vtkNewFoamIOobject io(this->CasePath,
      this->Parent->GetIsSinglePrecisionBinary() != 0);
  vtkStdString filePath(path + fileName);
  if (io.Open(filePath) || io.Open(filePath + ".gz"))
    {
    io.Close();
    // set points/faces location to current timesteps value
    polyMeshDir->SetValue(timeI, this->TimeNames->GetValue(timeI));
    }
  else
    {
    if (timeI != 0)
      {
      // set points/faces location to previous timesteps value
      polyMeshDir->SetValue(timeI, polyMeshDir->GetValue(timeI - 1));
      }
    else
      {
      filePath = this->CasePath + "constant" + this->RegionPath()
          + "/polyMesh/" + fileName;
      if(io.Open(filePath) || io.Open(filePath + ".gz"))
        {
        // set points/faces to constant
        polyMeshDir->SetValue(timeI, "constant");
        }
      }
    }
}

//-----------------------------------------------------------------------------
// create a Lookup Table containing the location of the points
// and faces files for each time steps mesh
void vtkNewOpenFOAMReaderPrivate::PopulatePolyMeshDirArrays()
{
  // intialize size to number of timesteps
  const int nSteps = this->TimeValues->GetNumberOfTuples();
  this->PolyMeshPointsDir->SetNumberOfValues(nSteps);
  this->PolyMeshFacesDir->SetNumberOfValues(nSteps);

  // loop through each timestep
  for (int i = 0; i < nSteps; i++)
    {
    // create the path to the timestep
    vtkStdString polyMeshPath = this->TimeRegionPath(i) + "/polyMesh/";
    AppendMeshDirToArray(this->PolyMeshPointsDir, polyMeshPath, "points", i);
    AppendMeshDirToArray(this->PolyMeshFacesDir, polyMeshPath, "faces", i);
    }
  return;
}

//-----------------------------------------------------------------------------
// read the points file into a vtkFloatArray
vtkFloatArray* vtkNewOpenFOAMReaderPrivate::ReadPointsFile()
{
  if (this->PolyMeshPointsDir->GetValue(this->TimeStep) == "")
    {
    vtkErrorMacro(<<"Cannot find path to points file for region \""
        << this->RegionName.c_str() << "\" in case " << this->CasePath.c_str());
    return NULL;
    }

  // path to points file
  const vtkStdString pointPath =
      this->CurrentTimeRegionMeshPath(this->PolyMeshPointsDir) + "points";

  vtkNewFoamIOobject io(this->CasePath,
      this->Parent->GetIsSinglePrecisionBinary() != 0);
  if (!(io.Open(pointPath) || io.Open(pointPath + ".gz")))
    {
    vtkErrorMacro(<<"Error opening " << io.GetFileName().c_str() << ": "
        << io.GetError().c_str());
    return NULL;
    }

  vtkNewFoamEntryValue dict(NULL);
  try
    {
    dict.ReadNonuniformList<vtkNewFoamToken::VECTORLIST,
    vtkNewFoamEntryValue::vectorListTraits<vtkFloatArray, float, 3, false> >(io);
    }
  catch(vtkNewFoamError& e)
    {
    vtkErrorMacro(<<"Error reading line " << io.GetLineNumber()
        << " of " << io.GetFileName().c_str() << ": " << e.c_str());
    return NULL;
    }

  vtkFloatArray *pointArray = static_cast<vtkFloatArray *>(dict.Ptr());

  // set the number of points
  this->NumPoints = pointArray->GetNumberOfTuples();

  return pointArray;
}

//-----------------------------------------------------------------------------
// read the faces into a vtkNewFoamIntVectorVector
vtkNewFoamIntVectorVector * vtkNewOpenFOAMReaderPrivate::ReadFacesFile(
    const vtkStdString &facePathIn)
{
  const vtkStdString facePath(facePathIn + "faces");

  vtkNewFoamIOobject io(this->CasePath,
      this->Parent->GetIsSinglePrecisionBinary() != 0);
  if (!(io.Open(facePath) || io.Open(facePath + ".gz")))
    {
    vtkErrorMacro(<<"Error opening " << io.GetFileName().c_str() << ": "
        << io.GetError().c_str() << ". If you are trying to read a parallel "
        "decomposed case, set Case Type to Decomposed Case.");
    return NULL;
    }

  vtkNewFoamEntryValue dict(NULL);
  try
    {
    dict.ReadLabelListList(io);
    }
  catch(vtkNewFoamError& e)
    {
    vtkErrorMacro(<<"Error reading line " << io.GetLineNumber()
        << " of " << io.GetFileName().c_str() << ": " << e.c_str());
    return NULL;
    }
  return static_cast<vtkNewFoamIntVectorVector *>(dict.Ptr());
}

//-----------------------------------------------------------------------------
// read the owner and neighbor file and create cellFaces
vtkNewFoamIntVectorVector * vtkNewOpenFOAMReaderPrivate::ReadOwnerNeighborFiles(
    const vtkStdString &ownerNeighborPath, vtkNewFoamIntVectorVector *facePoints)
{
  vtkNewFoamIOobject io(this->CasePath,
      this->Parent->GetIsSinglePrecisionBinary() != 0);
  vtkStdString ownerPath(ownerNeighborPath + "owner");
  if (io.Open(ownerPath) || io.Open(ownerPath + ".gz"))
    {
    vtkNewFoamEntryValue ownerDict(NULL);
    try
      {
      ownerDict.ReadNonuniformList<vtkNewFoamToken::LABELLIST,
      vtkNewFoamEntryValue::listTraits<vtkIntArray, int> >(io);
      }
    catch(vtkNewFoamError& e)
      {
      vtkErrorMacro(<<"Error reading line " << io.GetLineNumber()
          << " of " << io.GetFileName().c_str() << ": " << e.c_str());
      return NULL;
      }

    io.Close();

    const vtkStdString neighborPath(ownerNeighborPath + "neighbour");
    if (!(io.Open(neighborPath) || io.Open(neighborPath + ".gz")))
      {
      vtkErrorMacro(<<"Error opening " << io.GetFileName().c_str() << ": "
          << io.GetError().c_str());
      return NULL;
      }

    vtkNewFoamEntryValue neighborDict(NULL);
    try
      {
      neighborDict.ReadNonuniformList<vtkNewFoamToken::LABELLIST,
      vtkNewFoamEntryValue::listTraits<vtkIntArray, int> >(io);
      }
    catch(vtkNewFoamError& e)
      {
      vtkErrorMacro(<<"Error reading line " << io.GetLineNumber()
          << " of " << io.GetFileName().c_str() << ": " << e.c_str());
      return NULL;
      }

    this->FaceOwner = static_cast<vtkIntArray *>(ownerDict.Ptr());
    vtkIntArray &faceOwner = *this->FaceOwner;
    vtkIntArray &faceNeighbor = neighborDict.LabelList();

    const int nFaces = faceOwner.GetNumberOfTuples();
    const int nNeiFaces = faceNeighbor.GetNumberOfTuples();

    if (nFaces < nNeiFaces)
      {
      vtkErrorMacro(<<"Number of owner faces " << nFaces
          << " must be equal or larger than number of neighbor faces "
          << nNeiFaces);
      return NULL;
      }

    // the number of owner faces can be smaller than the number of
    // faces in sliding mesh cases
    if (nFaces > facePoints->GetNumberOfElements())
      {
      vtkWarningMacro(<<"Number of faces in faces "
          << facePoints->GetNumberOfElements()
          << " must be equal or larger than number of owner faces "
          << nFaces);
      return NULL;
      }

    // add the face numbers to the correct cell cf. Terry's code and
    // src/OpenFOAM/meshes/primitiveMesh/primitiveMeshCells.C
    // find the number of cells
    int nCells = -1;
    for (int faceI = 0; faceI < nNeiFaces; faceI++)
      {
      const int ownerCell = faceOwner.GetValue(faceI);
      if (nCells < ownerCell) // max(nCells, faceOwner[i])
        {
        nCells = ownerCell;
        }
      // we do need to take neighbor faces into account since all the
      // surrounding faces of a cell can be neighbors for a valid mesh
      const int neighborCell = faceNeighbor.GetValue(faceI);
      if (nCells < neighborCell) // max(nCells, faceNeighbor[i])
        {
        nCells = neighborCell;
        }
      }
    for (int faceI = nNeiFaces; faceI < nFaces; faceI++)
      {
      const int ownerCell = faceOwner.GetValue(faceI);
      if (nCells < ownerCell) // max(nCells, faceOwner[i])
        {
        nCells = ownerCell;
        }
      }
    nCells++;

    if (nCells == 0)
      {
      vtkWarningMacro(<<"The mesh contains no cells");
      }

    // set the number of cells
    this->NumCells = nCells;

    // create cellFaces with the length of the body undetermined
    vtkNewFoamIntVectorVector *cells = new vtkNewFoamIntVectorVector(nCells, 1);

    // count number of faces for each cell
    int *cfiPtr = cells->GetIndices()->GetPointer(0);
    for (int cellI = 0; cellI <= nCells; cellI++)
      {
      cfiPtr[cellI] = 0;
      }
    int nTotalCellFaces = 0;
    cfiPtr++; // offset +1
    for (int faceI = 0; faceI < nNeiFaces; faceI++)
      {
      const int ownerCell = faceOwner.GetValue(faceI);
      // simpleFoam/pitzDaily3Blocks has faces with owner cell number -1
      if (ownerCell >= 0)
        {
        cfiPtr[ownerCell]++;
        nTotalCellFaces++;
        }
      const int neighborCell=faceNeighbor.GetValue(faceI);
      if (neighborCell >= 0)
        {
        cfiPtr[neighborCell]++;
        nTotalCellFaces++;
        }
      }
    for (int faceI = nNeiFaces; faceI < nFaces; faceI++)
      {
      const int ownerCell = faceOwner.GetValue(faceI);
      if (ownerCell >= 0)
        {
        cfiPtr[ownerCell]++;
        nTotalCellFaces++;
        }
      }
    cfiPtr--; // revert offset +1

    // allocate cellFaces. To reduce the numbers of new/delete operations we
    // allocate memory space for all faces linearly
    cells->ResizeBody(nTotalCellFaces);

    // accumulate the number of cellFaces to create cellFaces indices
    // and copy them to a temporary array
    vtkIntArray *tmpFaceIndices = vtkIntArray::New();
    tmpFaceIndices->SetNumberOfValues(nCells + 1);
    int *tfiPtr = tmpFaceIndices->GetPointer(0);
    tfiPtr[0] = 0;
    for (int cellI = 1; cellI <= nCells; cellI++)
      {
      tfiPtr[cellI] = (cfiPtr[cellI] += cfiPtr[cellI - 1]);
      }

    // add face numbers to cell-faces list
    vtkIntArray *cellFacesList = cells->GetBody();
    for (int faceI = 0; faceI < nNeiFaces; faceI++)
      {
      const int ownerCell = faceOwner.GetValue(faceI); // must be a signed int
      // simpleFoam/pitzDaily3Blocks has faces with owner cell number -1
      if (ownerCell >= 0)
        {
        cellFacesList->SetValue(tfiPtr[ownerCell]++, faceI);
        }
      const int neighborCell = faceNeighbor.GetValue(faceI);
      if (neighborCell >= 0)
        {
        cellFacesList->SetValue(tfiPtr[neighborCell]++, faceI);
        }
      }
    for (int faceI = nNeiFaces; faceI < nFaces; faceI++)
      {
      const int ownerCell = faceOwner.GetValue(faceI); // must be a signed int
      // simpleFoam/pitzDaily3Blocks has faces with owner cell number -1
      if (ownerCell >= 0)
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
    if (!(io.Open(cellsPath) || io.Open(cellsPath + ".gz")))
      {
      vtkErrorMacro(<<"Error opening " << io.GetFileName().c_str() << ": "
          << io.GetError().c_str());
      return NULL;
      }
    vtkNewFoamEntryValue cellsDict(NULL);
    try
      {
      cellsDict.ReadLabelListList(io);
      }
    catch(vtkNewFoamError& e)
      {
      vtkErrorMacro(<<"Error reading line " << io.GetLineNumber()
          << " of " << io.GetFileName().c_str() << ": " << e.c_str());
      return NULL;
      }

    vtkNewFoamIntVectorVector *cells =
        static_cast<vtkNewFoamIntVectorVector *>(cellsDict.Ptr());
    this->NumCells = cells->GetNumberOfElements();
    const int nFaces = facePoints->GetNumberOfElements();

    // create face owner list
    this->FaceOwner = vtkIntArray::New();
    this->FaceOwner->SetNumberOfTuples(nFaces);
    for (int faceI = 0; faceI < nFaces; faceI++)
      {
      this->FaceOwner->SetValue(faceI, -1);
      }
    for (int cellI = 0; cellI < this->NumCells; cellI++)
      {
      const int nCellFaces = cells->GetSize(cellI);
      const int *cellFaces = cells->operator[](cellI);
      for (int faceI = 0; faceI < nCellFaces; faceI++)
        {
        const int f = cellFaces[faceI];
        if (f < 0 || f >= nFaces) // make sure the face number is valid
          {
          vtkErrorMacro("Face number " << f << " in cell " << cellI
              << " exceeds the number of faces " << nFaces);
          this->FaceOwner->Delete();
          this->FaceOwner = NULL;
          delete cells;
          return NULL;
          }
        const int owner = this->FaceOwner->GetValue(f);
        if (owner == -1 || owner > cellI)
          {
          this->FaceOwner->SetValue(f, cellI);
          }
        }
      }
    // check for unused faces
    for (int faceI = 0; faceI < nFaces; faceI++)
      {
      if (this->FaceOwner->GetValue(faceI) == -1)
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
bool vtkNewOpenFOAMReaderPrivate::CheckFacePoints(
    vtkNewFoamIntVectorVector *facePoints)
{
  const int nFaces = facePoints->GetNumberOfElements();

  for (int faceI = 0; faceI < nFaces; faceI++)
    {
    const int nPoints = facePoints->GetSize(faceI);
    const int *pointList = facePoints->operator[](faceI);
    if (nPoints < 3)
      {
      vtkErrorMacro(<< "Face " << faceI << " has only " << nPoints
          << " points which is not enough to constitute a face"
          " (a face must have at least 3 points)");
      return false;
      }
    for (int pointI = 0; pointI < nPoints; pointI++)
      {
      const int p = pointList[pointI];
      if (p < 0 || p >= this->NumPoints)
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
void vtkNewOpenFOAMReaderPrivate::InsertCellsToGrid(
    vtkUnstructuredGrid* internalMesh,
    const vtkNewFoamIntVectorVector *cellsFaces,
    const vtkNewFoamIntVectorVector *facesPoints, vtkFloatArray *pointArray,
    vtkIdTypeArray *additionalCells, vtkIntArray *cellList)
{
  const int maxNPoints = 128; // assume max number of points per cell
  vtkIdList* cellPoints = vtkIdList::New();
  cellPoints->SetNumberOfIds(maxNPoints);
  const int nCells = (cellList == NULL ? this->NumCells
      : cellList->GetNumberOfTuples());
  int nAdditionalPoints = 0;
  this->NumTotalAdditionalCells = 0;

  // alias
  const vtkNewFoamIntVectorVector& facePoints = *facesPoints;

  for (int cellI = 0; cellI < nCells; cellI++)
    {
    int cellId;
    if (cellList == NULL)
      {
      cellId = cellI;
      }
    else
      {
      cellId = cellList->GetValue(cellI);
      if (cellId >= this->NumCells)
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
    const int nCellFaces = cellsFaces->GetSize(cellId);

    // determine type of the cell
    // cf. src/OpenFOAM/meshes/meshShapes/cellMatcher/{hex|prism|pyr|tet}-
    // Matcher.C
    int cellType = VTK_CONVEX_POINT_SET;
    if (nCellFaces == 6)
      {
      int j = 0;
      for (; j < nCellFaces; j++)
        {
        if (facePoints.GetSize(cellFaces[j]) != 4)
          {
          break;
          }
        }
      if (j == nCellFaces)
        {
        cellType = VTK_HEXAHEDRON;
        }
      }
    else if (nCellFaces == 5)
      {
      int nTris = 0, nQuads = 0;
      for (int j = 0; j < nCellFaces; j++)
        {
        const int nPoints = facePoints.GetSize(cellFaces[j]);
        if (nPoints == 3)
          {
          nTris++;
          }
        else if (nPoints == 4)
          {
          nQuads++;
          }
        else
          {
          break;
          }
        }
      if (nTris == 2 && nQuads == 3)
        {
        cellType = VTK_WEDGE;
        }
      else if (nTris == 4 && nQuads == 1)
        {
        cellType = VTK_PYRAMID;
        }
      }
    else if (nCellFaces == 4)
      {
      int j = 0;
      for (; j < nCellFaces; j++)
        {
        if (facePoints.GetSize(cellFaces[j]) != 3)
          {
          break;
          }
        }
      if (j == nCellFaces)
        {
        cellType = VTK_TETRA;
        }
      }

    // not a Hex/Wedge/Pyramid/Tetra
    if (cellType == VTK_CONVEX_POINT_SET)
      {
      int nPoints = 0;
      for (int j = 0; j < nCellFaces; j++)
        {
        nPoints += facePoints.GetSize(cellFaces[j]);
        }
      if (nPoints == 0)
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

      if (this->FaceOwner->GetValue(cellBaseFaceId) == cellId)
        {
        // if it is an owner face flip the points
        for (int j = 0; j < 4; j++)
          {
          cellPoints->SetId(j, face0Points[3 - j]);
          }
        }
      else
        {
        // add base face to cell points
        for (int j = 0; j < 4; j++)
          {
          cellPoints->SetId(j, face0Points[j]);
          }
        }
      const int baseFacePoint0 = cellPoints->GetId(0);
      const int baseFacePoint2 = cellPoints->GetId(2);
      int cellOppositeFaceI = -1, pivotPoint = -1;
      int dupPoint = -1;
      for (int faceI = 1; faceI < 5; faceI++) // skip face 0 and 5
        {
        const int cellFaceI = cellFaces[faceI];
        const int *faceIPoints = facePoints[cellFaceI];
        int foundDup = -1, pointI = 0;
        for (; pointI < 4; pointI++) // each point
          {
          const int faceIPointI = faceIPoints[pointI];
          // matching two points in base face is enough to find a
          // duplicated point since neighboring faces share two
          // neighboring points (i. e. an edge)
          if (baseFacePoint0 == faceIPointI)
            {
            foundDup = 0;
            break;
            }
          else if (baseFacePoint2 == faceIPointI)
            {
            foundDup = 2;
            break;
            }
          }
        if (foundDup >= 0)
          {
          // find the pivot point if still haven't
          if (pivotPoint == -1)
            {
            dupPoint = foundDup;

            const int faceINextPoint = faceIPoints[(pointI + 1) % 4];

            // if the next point of the faceI-th face matches the
            // previous point of the base face use the previous point
            // of the faceI-th face as the pivot point; or use the
            // next point otherwise
            if (faceINextPoint == (this->FaceOwner->GetValue(cellFaceI)
                == cellId ? cellPoints->GetId(1 + foundDup)
                : cellPoints->GetId(3 - foundDup)))
              {
              pivotPoint = faceIPoints[(3 + pointI) % 4];
              }
            else
              {
              pivotPoint = faceINextPoint;
              }

            if (cellOppositeFaceI >= 0)
              {
              break;
              }
            }
          }
        else
          {
          // if no duplicated point found, faceI is the opposite face
          cellOppositeFaceI = cellFaceI;

          if (pivotPoint >= 0)
            {
            break;
            }
          }
        }

      // if the opposite face is not found until face 4, face 5 is
      // always the opposite face
      if (cellOppositeFaceI == -1)
        {
        cellOppositeFaceI = cellFaces[5];
        }

      // find the pivot point in opposite face
      const int *oppositeFacePoints = facePoints[cellOppositeFaceI];
      int pivotPointI = 0;
      for (; pivotPointI < 4; pivotPointI++)
        {
        if (oppositeFacePoints[pivotPointI] == pivotPoint)
          {
          break;
          }
        }

      // shift the pivot point if the point corresponds to point 2
      // of the base face
      if (dupPoint == 2)
        {
        pivotPointI = (pivotPointI + 2) % 4;
        }
      // copy the face-point list of the opposite face to cell-point list
      int basePointI = 4;
      if (this->FaceOwner->GetValue(cellOppositeFaceI) == cellId)
        {
        for (int pointI = pivotPointI; pointI < 4; pointI++)
          {
          cellPoints->SetId(basePointI++, oppositeFacePoints[pointI]);
          }
        for (int pointI = 0; pointI < pivotPointI; pointI++)
          {
          cellPoints->SetId(basePointI++, oppositeFacePoints[pointI]);
          }
        }
      else
        {
        for (int pointI = pivotPointI; pointI >= 0; pointI--)
          {
          cellPoints->SetId(basePointI++, oppositeFacePoints[pointI]);
          }
        for (int pointI = 3; pointI > pivotPointI; pointI--)
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
      for (int j = 0; j < 5; j++)
        {
        if (facePoints.GetSize(cellFaces[j]) == 3)
          {
          baseFaceId = j;
          break;
          }
        }

      // get first face in correct order
      const int cellBaseFaceId = cellFaces[baseFaceId];
      const int *face0Points = facePoints[cellBaseFaceId];

      if (this->FaceOwner->GetValue(cellBaseFaceId) == cellId)
        {
        for (int j = 0; j < 3; j++)
          {
          cellPoints->SetId(j, face0Points[j]);
          }
        }
      else
        {
        // if it is a neighbor face flip the points
        for (int j = 0; j < 3; j++)
          {
          // add base face to cell points
          cellPoints->SetId(j, face0Points[2 - j]);
          }
        }
      const int baseFacePoint0 = cellPoints->GetId(0);
      const int baseFacePoint2 = cellPoints->GetId(2);
      int cellOppositeFaceI = -1, pivotPoint = -1;
      bool dupPoint2 = false;
      for (int faceI = 0; faceI < 5; faceI++)
        {
        if (faceI == baseFaceId)
          {
          continue;
          }
        const int cellFaceI = cellFaces[faceI];
        if (facePoints.GetSize(cellFaceI) == 3)
          {
          cellOppositeFaceI = cellFaceI;
          }
        // find the pivot point if still haven't
        else if (pivotPoint == -1)
          {
          const int *faceIPoints = facePoints[cellFaceI];
          bool found0Dup = false, found2Dup = false;
          int pointI = 0;
          for (; pointI < 4; pointI++) // each point
            {
            const int faceIPointI = faceIPoints[pointI];
            // matching two points in base face is enough to find a
            // duplicated point since neighboring faces share two
            // neighboring points (i. e. an edge)
            if (baseFacePoint0 == faceIPointI)
              {
              found0Dup = true;
              break;
              }
            else if (baseFacePoint2 == faceIPointI)
              {
              found2Dup = true;
              break;
              }
            }
          // the matching point must always be found so omit the check
          int baseFacePrevPoint, baseFaceNextPoint;
          if (found0Dup)
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
          if (faceINextPoint
              == (this->FaceOwner->GetValue(cellFaceI) == cellId ? baseFacePrevPoint
                  : baseFaceNextPoint))
            {
            pivotPoint = faceIPrevPoint;
            }
          else
            {
            pivotPoint = faceINextPoint;
            }
          }

        // break when both of opposite face and pivot point are found
        if (cellOppositeFaceI >= 0 && pivotPoint >= 0)
          {
          break;
          }
        }

      // find the pivot point in opposite face
      const int *oppositeFacePoints = facePoints[cellOppositeFaceI];
      int pivotPointI = 0;
      for (; pivotPointI < 3; pivotPointI++)
        {
        if (oppositeFacePoints[pivotPointI] == pivotPoint)
          {
          break;
          }
        }

      if (this->FaceOwner->GetValue(cellOppositeFaceI) == cellId)
        {
        if (dupPoint2)
          {
          pivotPointI = (pivotPointI + 2) % 3;
          }
        int basePointI = 3;
        for (int pointI = pivotPointI; pointI >= 0; pointI--)
          {
          cellPoints->SetId(basePointI++, oppositeFacePoints[pointI]);
          }
        for (int pointI = 2; pointI > pivotPointI; pointI--)
          {
          cellPoints->SetId(basePointI++, oppositeFacePoints[pointI]);
          }
        }
      else
        {
        // shift the pivot point if the point corresponds to point 2
        // of the base face
        if (dupPoint2)
          {
          pivotPointI = (1 + pivotPointI) % 3;
          }
        // copy the face-point list of the opposite face to cell-point list
        int basePointI = 3;
        for (int pointI = pivotPointI; pointI < 3; pointI++)
          {
          cellPoints->SetId(basePointI++, oppositeFacePoints[pointI]);
          }
        for (int pointI = 0; pointI < pivotPointI; pointI++)
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
      if (cellType == VTK_PYRAMID)
        {
        for (int j = 0; j < nCellFaces; j++)
          {
          if (facePoints.GetSize(cellFaces[j]) == 4)
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
      const vtkIdType nBaseFacePoints = facePoints.GetSize(cellBaseFaceId);
      if (this->FaceOwner->GetValue(cellBaseFaceId) == cellId)
        {
        // if it is an owner face flip the points
        for (vtkIdType j = 0; j < nBaseFacePoints; j++)
          {
          cellPoints->SetId(j, baseFacePoints[nBaseFacePoints - 1 - j]);
          }
        }
      else
        {
        for (vtkIdType j = 0; j < nBaseFacePoints; j++)
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
      for (vtkIdType j = 0; j < nBaseFacePoints; j++)
        {
        // if point 1 of the adjacent face matches point j of the base face...
        if (cellPoints->GetId(j) == adjacentFacePoint1)
          {
          // if point 2 of the adjacent face matches the previous point
          // of the base face use point 0 of the adjacent face as the
          // pivot point; use point 2 otherwise
          cellPoints->SetId(
              nBaseFacePoints,
              (adjacentFacePoints[2]
                  == cellPoints->GetId((this->FaceOwner->GetValue(cellAdjacentFaceId)
                      == cellId ? (j + 1) : (nBaseFacePoints + j - 1))
                      % nBaseFacePoints)) ? adjacentFacePoints[0]
                  : adjacentFacePoints[2]);
          foundDup = true;
          break;
          }
        }
      // if point 1 of the adjacent face does not match any points of
      // the base face, it's the pivot point
      if (!foundDup)
        {
        cellPoints->SetId(nBaseFacePoints, adjacentFacePoint1);
        }

      // create the tetra cell and insert it into the mesh
      internalMesh->InsertNextCell(cellType, nPoints, cellPoints->GetPointer(0));
      }

    // erronous cells
    else if (cellType == VTK_EMPTY_CELL)
      {
      vtkWarningMacro("Warning: No points in cellId " << cellId);
      internalMesh->InsertNextCell(VTK_EMPTY_CELL, 0, cellPoints->GetPointer(0));
      }

    // OFpolyhedron || vtkConvexPointSet
    else
      {
      if (additionalCells != NULL) // decompose into tets and pyramids
        {
        // calculate cell centroid and insert it to point list
        this->AdditionalCellPoints->push_back(vtkIntArray::New());
        vtkIntArray *polyCellPoints = this->AdditionalCellPoints->back();
        float centroid[3];
        centroid[0] = centroid[1] = centroid[2] = 0.0F;
        for (int j = 0; j < nCellFaces; j++)
          {
          // remove duplicate points from faces
          const int cellFacesJ = cellFaces[j];
          const int *faceJPoints = facePoints[cellFacesJ];
          const size_t nFaceJPoints = facePoints.GetSize(cellFacesJ);
          for (size_t k = 0; k < nFaceJPoints; k++)
            {
            const int faceJPointK = faceJPoints[k];
            bool foundDup = false;
            for (vtkIdType l = 0; l < static_cast<vtkIdType>(polyCellPoints->GetDataSize()); l++)
              {
              if (polyCellPoints->GetValue(l) == faceJPointK)
                {
                foundDup = true;
                break; // look no more
                }
              }
            if (!foundDup)
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
        const float weight = 1.0F
            / static_cast<float>(polyCellPoints->GetDataSize());
        centroid[0] *= weight;
        centroid[1] *= weight;
        centroid[2] *= weight;
        pointArray->InsertNextTuple(centroid);

        // polyhedron decomposition.
        // a tweaked algorithm based on applications/utilities/postProcessing/
        // graphics/PVFoamReader/vtkNewFoam/vtkNewFoamAddInternalMesh.C
        bool insertDecomposedCell = true;
        int nAdditionalCells = 0;
        for (int j = 0; j < nCellFaces; j++)
          {
          const int cellFacesJ = cellFaces[j];
          const int *faceJPoints = facePoints[cellFacesJ];
          const int nFaceJPoints = facePoints.GetSize(cellFacesJ);
          const int flipNeighbor = (this->FaceOwner->GetValue(cellFacesJ)
              == cellId ? -1 : 1);
          const int nTris = nFaceJPoints % 2;

          int vertI = 2;

          // shift the start and end of the vertex loop if the
          // triangle of a decomposed face is going to be flat. Far
          // from perfect but better than nothing to avoid flat cells
          // which stops time integration of Stream Tracer especially
          // for split-hex unstructured meshes created by
          // e. g. autoRefineMesh
          if (nFaceJPoints >= 5 && nTris)
            {
            float *point0, *point1, *point2;
            point0 = pointArray->GetPointer(3 * faceJPoints[nFaceJPoints - 1]);
            point1 = pointArray->GetPointer(3 * faceJPoints[0]);
            point2 = pointArray->GetPointer(3 * faceJPoints[nFaceJPoints - 2]);
            float vsizeSqr1 = 0.0F, vsizeSqr2 = 0.0F, dotProduct = 0.0F;
            for (int i = 0; i < 3; i++)
              {
              const float v1 = point1[i] - point0[i], v2 = point2[i]
                  - point0[i];
              vsizeSqr1 += v1 * v1;
              vsizeSqr2 += v2 * v2;
              dotProduct += v1 * v2;
              }
            // compare in squared representation to avoid using sqrt()
            if (dotProduct * (float) fabs(dotProduct) / (vsizeSqr1 * vsizeSqr2)
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
          for (; vertI < nQuadVerts; vertI += 2)
            {
            cellPoints->SetId(1, faceJPoints[vertI - flipNeighbor]);
            cellPoints->SetId(2, faceJPoints[vertI]);
            cellPoints->SetId(3, faceJPoints[vertI + flipNeighbor]);

            // if the decomposed cell is the first one insert it to
            // the original position; or append to the decomposed cell
            // list otherwise
            if (insertDecomposedCell)
              {
              internalMesh->InsertNextCell(VTK_PYRAMID, 5,
                  cellPoints->GetPointer(0));
              insertDecomposedCell = false;
              }
            else
              {
              nAdditionalCells++;
              additionalCells->InsertNextTupleValue(cellPoints->GetPointer(0));
              }
            }

          // if the number of vertices is odd there's a triangle
          if (nTris)
            {
            if (flipNeighbor == -1)
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

            if (insertDecomposedCell)
              {
              internalMesh->InsertNextCell(VTK_TETRA, 4,
                  cellPoints->GetPointer(0));
              insertDecomposedCell = false;
              }
            else
              {
              // set the 5th vertex number to -1 to distinguish a tetra cell
              cellPoints->SetId(4, -1);
              nAdditionalCells++;
              additionalCells->InsertNextTupleValue(cellPoints->GetPointer(0));
              }
            }
          }
        nAdditionalPoints++;
        this->AdditionalCellIds->InsertNextValue(cellId);
        this->NumAdditionalCells->InsertNextValue(nAdditionalCells);
        this->NumTotalAdditionalCells += nAdditionalCells;
        }
      else // don't decompose; use VTK_CONVEX_PONIT_SET
        {
        // get first face
        const int cellFaces0 = cellFaces[0];
        const int *baseFacePoints = facePoints[cellFaces0];
        const int nBaseFacePoints = facePoints.GetSize(cellFaces0);
        int nPoints = nBaseFacePoints;
        if (nPoints > maxNPoints)
          {
          vtkErrorMacro(<< "Too large polyhedron at cellId = " << cellId);
          return;
          }
        if (this->FaceOwner->GetValue(cellFaces0) == cellId)
          {
          // if it is an owner face flip the points
          // not sure if flipping is necessary but do it anyway
          for (int j = 0; j < nBaseFacePoints; j++)
            {
            cellPoints->SetId(j, baseFacePoints[nBaseFacePoints - 1 - j]);
            }
          }
        else
          {
          // add first face to cell points
          for (int j = 0; j < nBaseFacePoints; j++)
            {
            cellPoints->SetId(j, baseFacePoints[j]);
            }
          }

        // loop through faces and create a list of all points
        // j = 1 skip baseFace
        for (int j = 1; j < nCellFaces; j++)
          {
          // remove duplicate points from faces
          const int cellFacesJ = cellFaces[j];
          const int *faceJPoints = facePoints[cellFacesJ];
          const size_t nFaceJPoints = facePoints.GetSize(cellFacesJ);
          for (size_t k = 0; k < nFaceJPoints; k++)
            {
            const int faceJPointK = faceJPoints[k];
            bool foundDup = false;
            for (int l = 0; l < nPoints; l++)
              {
              if (cellPoints->GetId(l) == faceJPointK)
                {
                foundDup = true;
                break; // look no more
                }
              }
            if (!foundDup)
              {
              if (nPoints >= maxNPoints)
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
void vtkNewOpenFOAMReaderPrivate::SetBlockName(vtkMultiBlockDataSet *blocks,
    unsigned int blockI, const char *name)
{
  blocks->GetMetaData(blockI)->Set(vtkCompositeDataSet::NAME(), name);
}

//-----------------------------------------------------------------------------
// derive cell types and create the internal mesh
vtkUnstructuredGrid *vtkNewOpenFOAMReaderPrivate::MakeInternalMesh(
    const vtkNewFoamIntVectorVector *cellsFaces,
    const vtkNewFoamIntVectorVector *facesPoints, vtkFloatArray *pointArray)
{
  // Create Mesh
  vtkUnstructuredGrid* internalMesh = vtkUnstructuredGrid::New();
  internalMesh->Allocate(this->NumCells);

  if (this->Parent->GetDecomposePolyhedra())
    {
    // for polyhedral decomposition
    this->AdditionalCellIds = vtkIntArray::New();
    this->NumAdditionalCells = vtkIntArray::New();
    this->AdditionalCellPoints = new vtkNewFoamIntArrayVector;

    vtkIdTypeArray *additionalCells = vtkIdTypeArray::New();
    additionalCells->SetNumberOfComponents(5); // accommodates tetra or pyramid

    this->InsertCellsToGrid(internalMesh, cellsFaces, facesPoints, pointArray,
        additionalCells, NULL);

    // for polyhedral decomposition
    pointArray->Squeeze();
    this->AdditionalCellIds->Squeeze();
    this->NumAdditionalCells->Squeeze();
    additionalCells->Squeeze();

    // insert decomposed cells into mesh
    const int nComponents = additionalCells->GetNumberOfComponents();
    const int nAdditionalCells = additionalCells->GetNumberOfTuples();
    for (int i = 0; i < nAdditionalCells; i++)
      {
      if (additionalCells->GetComponent(i, 4) == -1)
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
bool vtkNewOpenFOAMReaderPrivate::InsertFacesToGrid(vtkPolyData *boundaryMesh,
    const vtkNewFoamIntVectorVector *facesPoints, int startFace, int endFace,
    vtkIntArray *boundaryPointMap, vtkIdList *facePointsVtkId,
    vtkIntArray *labels, bool isLookupValue)
{
  vtkPolyData &bm = *boundaryMesh;

  for (int j = startFace; j < endFace; j++)
    {
    int faceId;
    if (labels == NULL)
      {
      faceId = j;
      }
    else
      {
      faceId = labels->GetValue(j);
      if (faceId >= facesPoints->GetNumberOfElements())
        {
        vtkErrorMacro(<<"faceLabels id " << faceId
            << " exceeds the number of faces "
            << facesPoints->GetNumberOfElements());
        return false;
        }
      }
    const int *facePoints = facesPoints->operator[](faceId);
    vtkIdType nFacePoints = facesPoints->GetSize(faceId);

    if (isLookupValue)
      {
      for (vtkIdType k = 0; k < nFacePoints; k++)
        {
        facePointsVtkId->SetId(k, boundaryPointMap->LookupValue(facePoints[k]));
        }
      }
    else
      {
      if (boundaryPointMap)
        {
        for (vtkIdType k = 0; k < nFacePoints; k++)
          {
          facePointsVtkId->SetId(k, boundaryPointMap->GetValue(facePoints[k]));
          }
        }
      else
        {
        for (vtkIdType k = 0; k < nFacePoints; k++)
          {
          facePointsVtkId->SetId(k, facePoints[k]);
          }
        }
      }

    // triangle
    if (nFacePoints == 3)
      {
      bm.InsertNextCell(VTK_TRIANGLE, 3, facePointsVtkId->GetPointer(0));
      }
    // quad
    else if (nFacePoints == 4)
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
  return true;
}

//-----------------------------------------------------------------------------
// returns requested boundary meshes
vtkMultiBlockDataSet *vtkNewOpenFOAMReaderPrivate::MakeBoundaryMesh(
    const vtkNewFoamIntVectorVector *facesPoints, vtkFloatArray* pointArray)
{
  const vtkIdType nBoundaries = static_cast<vtkIdType>(this->BoundaryDict.size());

  // do a consistency check of BoundaryDict
  int previousEndFace = -1;
  for (int boundaryI = 0; boundaryI < nBoundaries; boundaryI++)
    {
    const vtkNewFoamBoundaryEntry &beI = this->BoundaryDict[boundaryI];
    const int startFace = beI.StartFace;
    const int nFaces = beI.NFaces;
    if (startFace < 0 || nFaces < 0)
      {
      vtkErrorMacro(<<"Neither of startFace " << startFace << " nor nFaces "
          << nFaces << " can be nagative for patch " << beI.BoundaryName.c_str());
      return NULL;
      }
    if (previousEndFace >= 0 && previousEndFace != startFace)
      {
      vtkErrorMacro(<<"The end face number " << previousEndFace - 1
          << " of patch "
          << this->BoundaryDict[boundaryI - 1].BoundaryName.c_str()
          << " is not consistent with the start face number " << startFace
          << " of patch " << beI.BoundaryName.c_str());
      return NULL;
      }
    previousEndFace = startFace + nFaces;
    }
  if (previousEndFace > facesPoints->GetNumberOfElements())
    {
    vtkErrorMacro(<<"The end face number " << previousEndFace - 1
        << " of the last patch "
        << this->BoundaryDict[nBoundaries - 1].BoundaryName.c_str()
        << " exceeds the number of faces " << facesPoints->GetNumberOfElements());
    return NULL;
    }

  vtkMultiBlockDataSet *boundaryMesh = vtkMultiBlockDataSet::New();

  if (this->Parent->GetCreateCellToPoint())
    {
    const int boundaryStartFace =
        (this->BoundaryDict.size() > 0 ? this->BoundaryDict[0].StartFace : 0);
    this->AllBoundaries = vtkPolyData::New();
    this->AllBoundaries->Allocate(facesPoints->GetNumberOfElements()
        - boundaryStartFace);
    }
  this->BoundaryPointMap = new vtkNewFoamIntArrayVector;

  vtkIntArray *nBoundaryPointsList = vtkIntArray::New();
  nBoundaryPointsList->SetNumberOfValues(nBoundaries);

  // count the max number of points per face and the number of points
  // (with duplicates) in mesh
  int maxNFacePoints = 0;
  for (int boundaryI = 0; boundaryI < nBoundaries; boundaryI++)
    {
    const int startFace = this->BoundaryDict[boundaryI].StartFace;
    const int endFace = startFace + this->BoundaryDict[boundaryI].NFaces;
    int nPoints = 0;
    for (int j = startFace; j < endFace; j++)
      {
      const int nFacePoints = facesPoints->GetSize(j);
      nPoints += nFacePoints;
      if (nFacePoints > maxNFacePoints)
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
  if (this->Parent->GetCreateCellToPoint())
    {
    this->InternalPoints = vtkIntArray::New();
    this->InternalPoints->SetNumberOfValues(this->NumPoints);
    for (int pointI = 0; pointI < this->NumPoints; pointI++)
      {
      this->InternalPoints->SetValue(pointI, -1);
      }

    // mark boundary points as 0
    for (int boundaryI = 0; boundaryI < nBoundaries; boundaryI++)
      {
      const vtkNewFoamBoundaryEntry &beI = this->BoundaryDict[boundaryI];
      if (beI.BoundaryType == vtkNewFoamBoundaryEntry::PHYSICAL
          || beI.BoundaryType == vtkNewFoamBoundaryEntry::PROCESSOR)
        {
        const int startFace = beI.StartFace;
        const int endFace = startFace + beI.NFaces;

        for (int j = startFace; j < endFace; j++)
          {
          const int *facePoints = facesPoints->operator[](j);
          const int nFacePoints = facesPoints->GetSize(j);
          for (int k = 0; k < nFacePoints; k++)
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

  if (this->Parent->GetCreateCellToPoint())
    {
    // create global to AllBoundaries point map
    for (int pointI = 0; pointI < this->NumPoints; pointI++)
      {
      if (this->InternalPoints->GetValue(pointI) == 0)
        {
        this->InternalPoints->SetValue(pointI, nAllBoundaryPoints);
        nAllBoundaryPoints++;
        }
      }

    if (this->ProcessorName != "")
      {
      // initialize physical-processor boundary shared point list
      procCellList.resize(nAllBoundaryPoints);
      pointTypes = vtkIntArray::New();
      pointTypes->SetNumberOfTuples(nAllBoundaryPoints);
      for (int pointI = 0; pointI < nAllBoundaryPoints; pointI++)
        {
        pointTypes->SetValue(pointI, 0);
        }
      }
    }

  for (int boundaryI = 0; boundaryI < nBoundaries; boundaryI++)
    {
    const vtkNewFoamBoundaryEntry &beI = this->BoundaryDict[boundaryI];
    const int nFaces = beI.NFaces;
    const int startFace = beI.StartFace;
    const int endFace = startFace + nFaces;

    if (this->Parent->GetCreateCellToPoint() && (beI.BoundaryType
        == vtkNewFoamBoundaryEntry::PHYSICAL || beI.BoundaryType
        == vtkNewFoamBoundaryEntry::PROCESSOR ))
      {
      // add faces to AllBoundaries
      this->InsertFacesToGrid(this->AllBoundaries, facesPoints, startFace,
          endFace, this->InternalPoints, facePointsVtkId, NULL, false);

      if (this->ProcessorName != "")
        {
        // mark belonging boundary types and, if PROCESSOR, cell numbers
        const int abStartFace = beI.AllBoundariesStartFace;
        const int abEndFace = abStartFace + beI.NFaces;
        for (int faceI = abStartFace; faceI < abEndFace; faceI++)
          {
          vtkIdType nPoints;
          vtkIdType *points;
          this->AllBoundaries->GetCellPoints(faceI, nPoints, points);
          if (beI.BoundaryType == vtkNewFoamBoundaryEntry::PHYSICAL)
            {
            for (int pointI = 0; pointI < nPoints; pointI++)
              {
              *pointTypes->GetPointer(points[pointI])
                  |= vtkNewFoamBoundaryEntry::PHYSICAL;
              }
            }
          else // PROCESSOR
            {
            for (int pointI = 0; pointI < nPoints; pointI++)
              {
              const int pointJ = points[pointI];
              *pointTypes->GetPointer(pointJ)
                  |= vtkNewFoamBoundaryEntry::PROCESSOR;
              procCellList[pointJ].push_back(faceI);
              }
            }
          }
        }
      }

    // skip below if inactive
    if (!beI.IsActive)
      {
      continue;
      }

    // create the mesh
    const unsigned int activeBoundaryI = boundaryMesh->GetNumberOfBlocks();
    vtkPolyData *bm = vtkPolyData::New();
    boundaryMesh->SetBlock(activeBoundaryI, bm);

    // set the name of boundary
    this->SetBlockName(boundaryMesh, activeBoundaryI, beI.BoundaryName.c_str());

    bm->Allocate(nFaces);
    const int nBoundaryPoints = nBoundaryPointsList->GetValue(boundaryI);

    // create global to boundary-local point map and boundary points
    vtkIntArray *boundaryPointList = vtkIntArray::New();
    boundaryPointList->SetNumberOfValues(nBoundaryPoints);
    int pointI = 0;
    for (int j = startFace; j < endFace; j++)
      {
      const int *facePoints = facesPoints->operator[](j);
      int nFacePoints = facesPoints->GetSize(j);
      for (int k = 0; k < nFacePoints; k++)
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
    for (int j = 0; j < nBoundaryPoints; j++)
      {
      const int pointJ = boundaryPointList->GetValue(j);
      if (pointJ != oldPointJ)
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

  if (this->Parent->GetCreateCellToPoint())
    {
    this->AllBoundaries->Squeeze();
    this->AllBoundariesPointMap = vtkIntArray::New();
    vtkIntArray &abpMap = *this->AllBoundariesPointMap;
    abpMap.SetNumberOfValues(nAllBoundaryPoints);

    // create lists of internal points and AllBoundaries points
    int nInternalPoints = 0;
    for (int pointI = 0, allBoundaryPointI = 0; pointI < this->NumPoints; pointI++)
      {
      const int globalPointId = this->InternalPoints->GetValue(pointI);
      if (globalPointId == -1)
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
    if (nInternalPoints > 0)
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

    if (this->ProcessorName != "")
      {
      // remove links to processor boundary faces from point-to-cell
      // links of physical-processor shared points to avoid cracky seams
      // on fixedValue-type boundaries which are noticeable when all the
      // decomposed meshes are appended
      this->AllBoundaries->BuildLinks();
      for (int pointI = 0; pointI < nAllBoundaryPoints; pointI++)
        {
        if (pointTypes->GetValue(pointI) == (vtkNewFoamBoundaryEntry::PHYSICAL
            | vtkNewFoamBoundaryEntry::PROCESSOR))
          {
          const vtkstd::vector<int> &procCells = procCellList[pointI];
          for (size_t cellI = 0; cellI < procCellList[pointI].size(); cellI++)
            {
            this->AllBoundaries->RemoveReferenceToCell(pointI, procCells[cellI]);
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
void vtkNewOpenFOAMReaderPrivate::TruncateFaceOwner()
{
  const int boundaryStartFace =
      this->BoundaryDict.size() > 0 ? this->BoundaryDict[0].StartFace
          : this->FaceOwner->GetNumberOfTuples();
  // all the boundary faces
  const int nBoundaryFaces = this->FaceOwner->GetNumberOfTuples()
      - boundaryStartFace;
  memmove(this->FaceOwner->GetPointer(0),
      this->FaceOwner->GetPointer(boundaryStartFace), sizeof(int)
          * nBoundaryFaces);
  this->FaceOwner->Resize(nBoundaryFaces);
}

//-----------------------------------------------------------------------------
// this is necessary due to the strange vtkDataArrayTemplate::Resize()
// implementation when the array size is to be extended
template <typename T1, typename T2> bool vtkNewOpenFOAMReaderPrivate::ExtendArray(
    T1 *array, const int nTuples)
{
  const int newSize = nTuples * array->GetNumberOfComponents();
  void *ptr = malloc(newSize * array->GetDataTypeSize());
  if (ptr == NULL)
    {
    return false;
    }
  memmove(ptr, array->GetVoidPointer(0), array->GetDataSize()
      * array->GetDataTypeSize());
  array->SetArray(static_cast<T2 *>(ptr), newSize, 0);
  return true;
}

#if 0
//-----------------------------------------------------------------------------
void vtkNewOpenFOAMReaderPrivate::CalculateReciprocalDelta(
  const vtkNewFoamIntVectorVector *facesPoints)
{
  const int nBoundaries = this->BoundaryDict.size();
  vtkFloatArray *allPoints
    = vtkFloatArray::SafeDownCast(this->InternalMesh->GetPoints()->GetData());

  this->ReciprocalDelta = new vtkNewFoamFloatArrayVector;
  for(int boundaryI = 0; boundaryI < nBoundaries; boundaryI++)
    {
    const vtkNewFoamBoundaryEntry &beI = this->BoundaryDict[boundaryI];
    if(beI.BoundaryType != vtkNewFoamBoundaryEntry::PHYSICAL)
      {
      this->ReciprocalDelta->push_back(NULL);
      continue;
      }

    this->ReciprocalDelta->push_back(vtkFloatArray::New());
    vtkFloatArray *rd = this->ReciprocalDelta->back();
    rd->SetNumberOfTuples(beI.NFaces);
    const int boundaryStartFace = this->BoundaryDict[0].StartFace;
    const int startFace = beI.StartFace;
    const int endFace = startFace + beI.NFaces;
    for(int faceI = startFace; faceI < endFace; faceI++)
      {
      // calculate patchInternal cell centroid
      const int cellId = this->FaceOwner->GetValue(faceI - boundaryStartFace);
      const int polyCellId = this->AdditionalCellIds->LookupValue(cellId);
      float cellCentroid0, cellCentroid1, cellCentroid2;
      if(this->Parent->GetDecomposePolyhedra()
        && this->AdditionalCellPoints != NULL && polyCellId >= 0)
        {
        const float *cellCentroidPtr
          = allPoints->GetPointer(3 * (this->NumPoints + polyCellId));
        cellCentroid0 = cellCentroidPtr[0];
        cellCentroid1 = cellCentroidPtr[1];
        cellCentroid2 = cellCentroidPtr[2];
        }
      else
        {
        cellCentroid0 = cellCentroid1 = cellCentroid2 = 0.0F;
        vtkIdType nPoints, *pointIds;
        this->InternalMesh->GetCellPoints(cellId, nPoints, pointIds);
        for(int pointI = 0; pointI < nPoints; pointI++)
          {
          const float *point = allPoints->GetPointer(3 * pointIds[pointI]);
          cellCentroid0 += point[0];
          cellCentroid1 += point[1];
          cellCentroid2 += point[2];
          }
        const float weight = 1.0F / static_cast<float>(nPoints);
        cellCentroid0 *= weight;
        cellCentroid1 *= weight;
        cellCentroid2 *= weight;
        }

      // calculate face centroid
      float faceCentroid0 = 0.0F, faceCentroid1 = 0.0F, faceCentroid2 = 0.0F;
      const int *facePoints = facesPoints->operator[](faceI);
      const int nFacePoints = facesPoints->GetSize(faceI);
      for(int pointI = 0; pointI < nFacePoints; pointI++)
        {
        const float *point = allPoints->GetPointer(3 * facePoints[pointI]);
        faceCentroid0 += point[0];
        faceCentroid1 += point[1];
        faceCentroid2 += point[2];
        }
      const float weight = 1.0F / static_cast<float>(nFacePoints);
      faceCentroid0 *= weight;
      faceCentroid1 *= weight;
      faceCentroid2 *= weight;

      // calculate face normal
      const float *point0 = allPoints->GetPointer(3 * facePoints[0]);
      const float *point1 = allPoints->GetPointer(3 * facePoints[1]);
      const float v00 = point0[0] - point1[0];
      const float v01 = point0[1] - point1[1];
      const float v02 = point0[2] - point1[2];
      const float *point2 = allPoints->GetPointer(3 * facePoints[2]);
      const float v10 = point2[0] - point1[0];
      const float v11 = point2[1] - point1[1];
      const float v12 = point2[2] - point1[2];
      float n0 = v11 * v02 - v12 * v01;
      float n1 = v12 * v00 - v10 * v02;
      float n2 = v10 * v01 - v11 * v00;
      float nWeight = sqrt(n0 * n0 + n1 * n1 + n2 * n2);

      const float delta0 = faceCentroid0 - cellCentroid0;
      const float delta1 = faceCentroid1 - cellCentroid1;
      const float delta2 = faceCentroid2 - cellCentroid2;
      float delta;
      if(nWeight == 0.0F)
        {
        // use raw distance
        delta = sqrt(delta0 * delta0 + delta1 * delta1 + delta2 * delta2);
        }
      else
        {
        // use inner product
        delta = (delta0 * n0 + delta1 * n1 + delta2 * n2) / nWeight;
        }
      rd->SetValue(faceI - beI.StartFace, delta == 0.0F ? delta : 1.0F / delta);
      }
    }
  this->AdditionalCellIds->ClearLookup();
}
#endif

//-----------------------------------------------------------------------------
// move polyhedral cell centroids
vtkPoints *vtkNewOpenFOAMReaderPrivate::MoveInternalMesh(
    vtkUnstructuredGrid *internalMesh, vtkFloatArray *pointArray)
{
  if (this->Parent->GetDecomposePolyhedra())
    {
    const vtkIdType nAdditionalCells = static_cast<vtkIdType>(this->AdditionalCellPoints->size());
    this->ExtendArray<vtkFloatArray, float>(pointArray, this->NumPoints
        + nAdditionalCells);
    for (int i = 0; i < nAdditionalCells; i++)
      {
      vtkIntArray *polyCellPoints = this->AdditionalCellPoints->operator[](i);
      float centroid[3];
      centroid[0] = centroid[1] = centroid[2] = 0.0F;
      const int nCellPoints = polyCellPoints->GetDataSize();
      for (int j = 0; j < nCellPoints; j++)
        {
        float *pointK = pointArray->GetPointer(3 * polyCellPoints->GetValue(j));
        centroid[0] += pointK[0];
        centroid[1] += pointK[1];
        centroid[2] += pointK[2];
        }
      const float weight = (nCellPoints ? 1.0F
          / static_cast<float>(nCellPoints) : 0.0F);
      centroid[0] *= weight;
      centroid[1] *= weight;
      centroid[2] *= weight;
      pointArray->InsertTuple(this->NumPoints + i, centroid);
      }
    }
  if (internalMesh->GetPoints()->GetNumberOfPoints() != pointArray->GetNumberOfTuples())
    {
    vtkErrorMacro(<< "The numbers of points for old points "
        << internalMesh->GetPoints()->GetNumberOfPoints() << " and new points"
        << pointArray->GetNumberOfTuples() << " don't match");
    return NULL;
    }

  // instantiate the points class
  vtkPoints *points = vtkPoints::New();
  points->SetData(pointArray);
  internalMesh->SetPoints(points);
  return points;
}

//-----------------------------------------------------------------------------
// move boundary points
void vtkNewOpenFOAMReaderPrivate::MoveBoundaryMesh(
    vtkMultiBlockDataSet *boundaryMesh, vtkFloatArray *pointArray)
{
  for (vtkIdType boundaryI = 0, activeBoundaryI = 0; boundaryI
    < static_cast<vtkIdType>(this->BoundaryDict.size()); boundaryI++)
    {
    if (this->BoundaryDict[boundaryI].IsActive)
      {
      vtkIntArray *bpMap = this->BoundaryPointMap->operator[](activeBoundaryI);
      const int nBoundaryPoints = bpMap->GetNumberOfTuples();
      vtkFloatArray *boundaryPointArray = vtkFloatArray::New();
      boundaryPointArray->SetNumberOfComponents(3);
      boundaryPointArray->SetNumberOfTuples(nBoundaryPoints);
      for (int pointI = 0; pointI < nBoundaryPoints; pointI++)
        {
        boundaryPointArray->SetTuple(pointI, bpMap->GetValue(pointI),
            pointArray);
        }
      vtkPoints *boundaryPoints = vtkPoints::New();
      boundaryPoints->SetData(boundaryPointArray);
      boundaryPointArray->Delete();
      vtkPolyData::SafeDownCast(boundaryMesh->GetBlock(activeBoundaryI))
      ->SetPoints(boundaryPoints);
      boundaryPoints->Delete();
      activeBoundaryI++;
      }
    }
}

//-----------------------------------------------------------------------------
// as of now the function does not do interpolation, but do just averaging.
void vtkNewOpenFOAMReaderPrivate::InterpolateCellToPoint(vtkFloatArray *pData,
    vtkFloatArray *iData, vtkPointSet *mesh, vtkIntArray *pointList,
    const int nPoints)
{
  if (nPoints == 0)
    {
    return;
    }

  // a dummy call to let GetPointCells() build the cell links if still not built
  // (not using BuildLinks() since it always rebuild links)
  vtkIdList *pointCells = vtkIdList::New();
  mesh->GetPointCells(0, pointCells);
  pointCells->Delete();

  // since vtkPolyData and vtkUnstructuredGrid do not share common
  // overloaded GetCellLink() or GetPointCells() functions we have to
  // do a tedious task
  vtkUnstructuredGrid *ug = vtkUnstructuredGrid::SafeDownCast(mesh);
  vtkPolyData *pd = vtkPolyData::SafeDownCast(mesh);
  vtkCellLinks *cl = NULL;
  if (ug)
    {
    cl = ug->GetCellLinks();
    }

  const int nComponents = iData->GetNumberOfComponents();

  if (nComponents == 1)
    {
    // a special case with the innermost componentI loop unrolled
    float *tuples = iData->GetPointer(0);
    for (int pointI = 0; pointI < nPoints; pointI++)
      {
      const int pI = (pointList ? pointList->GetValue(pointI) : pointI);
      unsigned short nCells;
      vtkIdType *cells;
      if (cl)
        {
        const vtkCellLinks::Link &l = cl->GetLink(pI);
        nCells = l.ncells;
        cells = l.cells;
        }
      else
        {
        pd->GetPointCells(pI, nCells, cells);
        }
      // use double intermediate variable for precision
      double interpolatedValue = 0.0;
      for (int cellI = 0; cellI < nCells; cellI++)
        {
        interpolatedValue += tuples[cells[cellI]];
        }
      interpolatedValue = (nCells ? interpolatedValue
          / static_cast<double>(nCells) : 0.0);
      pData->SetValue(pI, interpolatedValue);
      }
    }
  else if (nComponents == 3)
    {
    // a special case with the innermost componentI loop unrolled
    float *pDataPtr = pData->GetPointer(0);
    for (int pointI = 0; pointI < nPoints; pointI++)
      {
      const int pI = (pointList ? pointList->GetValue(pointI) : pointI);
      unsigned short nCells;
      vtkIdType *cells;
      if (cl)
        {
        const vtkCellLinks::Link &l = cl->GetLink(pI);
        nCells = l.ncells;
        cells = l.cells;
        }
      else
        {
        pd->GetPointCells(pI, nCells, cells);
        }
      // use double intermediate variables for precision
      const double weight = (nCells ? 1.0 / static_cast<double>(nCells) : 0.0);
      double summedValue0 = 0.0, summedValue1 = 0.0, summedValue2 = 0.0;

      // hand unrolling
      for (int cellI = 0; cellI < nCells; cellI++)
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
    for (int pointI = 0; pointI < nPoints; pointI++)
      {
      const int pI = (pointList ? pointList->GetValue(pointI) : pointI);
      unsigned short nCells;
      vtkIdType *cells;
      if (cl)
        {
        const vtkCellLinks::Link &l = cl->GetLink(pI);
        nCells = l.ncells;
        cells = l.cells;
        }
      else
        {
        pd->GetPointCells(pI, nCells, cells);
        }
      // use double intermediate variables for precision
      const double weight = (nCells ? 1.0 / static_cast<double>(nCells) : 0.0);
      float *interpolatedValue = &pDataPtr[nComponents * pI];
      // a bit strange loop order but this works fastest
      for (int componentI = 0; componentI < nComponents; componentI++)
        {
        const float *tuple = iData->GetPointer(componentI);
        double summedValue = 0.0;
        for (int cellI = 0; cellI < nCells; cellI++)
          {
          summedValue += tuple[nComponents * cells[cellI]];
          }
        interpolatedValue[componentI] = weight * summedValue;
        }
      }
    }
}

//-----------------------------------------------------------------------------
bool vtkNewOpenFOAMReaderPrivate::ReadFieldFile(vtkNewFoamIOobject *ioPtr,
    vtkNewFoamDict *dictPtr, const vtkStdString &varName,
    vtkDataArraySelection *selection)
{
  const vtkStdString varPath(this->CurrentTimeRegionPath() + "/" + varName);

  // open the file
  vtkNewFoamIOobject &io = *ioPtr;
  if (!io.Open(varPath))
    {
    vtkErrorMacro(<<"Error opening " << io.GetFileName().c_str() << ": "
        << io.GetError().c_str());
    return false;
    }

  // if the variable is disabled on selection panel then skip it
  if (selection->ArrayExists(io.GetObjectName().c_str()) && !selection->ArrayIsEnabled(io.GetObjectName().c_str()))
    {
    return false;
    }

  // read the field file into dictionary
  vtkNewFoamDict &dict = *dictPtr;
  if (!dict.Read(io))
    {
    vtkErrorMacro(<<"Error reading line " << io.GetLineNumber()
        << " of " << io.GetFileName().c_str() << ": " << io.GetError().c_str());
    return false;
    }

  if (dict.GetType() != vtkNewFoamToken::DICTIONARY)
    {
    vtkErrorMacro(<<"File " << io.GetFileName().c_str()
        << "is not valid as a field file");
    return false;
    }
  return true;
}

//-----------------------------------------------------------------------------
vtkFloatArray *vtkNewOpenFOAMReaderPrivate::FillField(vtkNewFoamEntry *entryPtr,
    vtkNewFoamEntry *refEntryPtr, int nElements, vtkNewFoamIOobject *ioPtr,
    const vtkStdString &fieldType, const bool isUniformFixedValueBC)
    // the isUniformFixedValueBC argument can be determined in compile-time
    // thus better use a template, which unfortunately is not possible due to
    // broken gcc-3.3
{
  vtkFloatArray *data;
  vtkNewFoamEntry &entry = *entryPtr;
  const vtkStdString &className = ioPtr->GetClassName();

  vtkNewFoamToken::tokenType refEntryType;
  int nRefComponents;
  float refTuple[9];
  if (refEntryPtr != NULL)
    {
    if (refEntryPtr->FirstValue().GetType() == vtkNewFoamToken::SCALAR
        || refEntryPtr->FirstValue().GetType() == vtkNewFoamToken::LABEL)
      {
      refEntryType = vtkNewFoamToken::SCALAR;
      nRefComponents = 1;
      refTuple[0] = refEntryPtr->ToFloat();
      }
    else if (refEntryPtr->FirstValue().GetType() == vtkNewFoamToken::LABELLIST)
      {
      refEntryType = vtkNewFoamToken::SCALARLIST;
      vtkIntArray &ll = refEntryPtr->LabelList();
      nRefComponents = ll.GetNumberOfTuples();
      for (int componentI = 0; componentI < nRefComponents; componentI++)
        {
        refTuple[componentI] = static_cast<float>(ll.GetValue(componentI));
        }
      }
    else if (refEntryPtr->FirstValue().GetType() == vtkNewFoamToken::SCALARLIST)
      {
      refEntryType = vtkNewFoamToken::SCALARLIST;
      vtkFloatArray &sl = refEntryPtr->ScalarList();
      nRefComponents = sl.GetNumberOfTuples();
      for (int componentI = 0; componentI < nRefComponents; componentI++)
        {
        refTuple[componentI] = sl.GetValue(componentI);
        }
      }
    else
      {
      vtkErrorMacro(<< "Wrong referenceLevel type");
      return NULL;
      }
    }

  // the isUniformFixedValueBC argument is for uniformFixedValue B.C.
  if (entry.FirstValue().GetIsUniform() == vtkNewFoamEntryValue::UNIFORM
      || isUniformFixedValueBC)
    {
    if (entry.FirstValue().GetType() == vtkNewFoamToken::SCALAR || entry.FirstValue().GetType() == vtkNewFoamToken::LABEL)
      {
      float num = entry.ToFloat();
      if (refEntryPtr != NULL)
        {
        if (refEntryType == vtkNewFoamToken::SCALAR)
          {
          num += refTuple[0];
          }
        else
          {
          vtkErrorMacro(<<"Wrong referenceLevel type for uniform scalar field");
          return NULL;
          }
        }
      data = vtkFloatArray::New();
      data->SetNumberOfValues(nElements);
      for (int i = 0; i < nElements; i++)
        {
        data->SetValue(i, num);
        }
      }
    else
      {
      float tuple[9];
      int nComponents;
      // have to determine the type of vector
      if (entry.FirstValue().GetType() == vtkNewFoamToken::LABELLIST)
        {
        vtkIntArray &ll = entry.LabelList();
        nComponents = ll.GetNumberOfTuples();
        for (int componentI = 0; componentI < nComponents; componentI++)
          {
          tuple[componentI] = static_cast<float>(ll.GetValue(componentI));
          }
        }
      else if (entry.FirstValue().GetType() == vtkNewFoamToken::SCALARLIST)
        {
        vtkFloatArray& sl = entry.ScalarList();
        nComponents = sl.GetSize();
        for (int componentI = 0; componentI < nComponents; componentI++)
          {
          tuple[componentI] = sl.GetValue(componentI);
          }
        }
      else
        {
        vtkErrorMacro(<< "Wrong list type for uniform field");
        return NULL;
        }

      if (refEntryPtr != NULL)
        {
        if (refEntryType == vtkNewFoamToken::SCALARLIST && nRefComponents == nComponents)
          {
          for (int componentI = 0; componentI < nComponents; componentI++)
            {
            tuple[componentI] += refTuple[componentI];
            }
          }
        else
          {
          vtkErrorMacro(<< "Wrong referenceLevel type");
          return NULL;
          }
        }

      if ((fieldType == "SphericalTensorField" && nComponents == 1)
          || (fieldType == "VectorField" && nComponents == 3) || (fieldType
          == "SymmTensorField" && nComponents == 6) || (fieldType
          == "TensorField" && nComponents == 9))
        {
        data = vtkFloatArray::New();
        data->SetNumberOfComponents(nComponents);
        data->SetNumberOfTuples(nElements);
#if vtksys_DATE_STAMP_FULL >= 20080620
        // swap the components of symmTensor to match the component
        // names in paraview
        if (nComponents == 6)
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
        for (int i = 0; i < nElements; i++)
          {
          data->SetTuple(i, tuple);
          }
        }
      else
        {
        vtkErrorMacro(<< "Number of components and field class doesn't match "
                      << "for " << ioPtr->GetFileName().c_str() << ". class = " << className.c_str()
                      << ", nComponents = " << nComponents);
        return NULL;
        }
      }
    }
  else if (entry.FirstValue().GetIsUniform() == vtkNewFoamEntryValue::NONUNIFORM)
    {
    if ((fieldType == "ScalarField" && entry.FirstValue().GetType() == vtkNewFoamToken::SCALARLIST) || ((fieldType
        == "VectorField" || fieldType == "SphericalTensorField" || fieldType
        == "SymmTensorField" || fieldType == "TensorField")
        && entry.FirstValue().GetType() == vtkNewFoamToken::VECTORLIST))
      {
      const int nTuples = entry.ScalarList().GetNumberOfTuples();
      if (nTuples != nElements)
        {
        vtkErrorMacro(<<"Number of cells/points in mesh and field don't match: "
            << "mesh = " << nElements << ", field = " << nTuples);
        return NULL;
        }
      data = static_cast<vtkFloatArray *>(entry.Ptr());

      if (refEntryPtr != NULL)
        {
        const int nComponents = data->GetNumberOfComponents();
        if ((refEntryType == vtkNewFoamToken::SCALAR
            && entry.FirstValue().GetType() == vtkNewFoamToken::SCALARLIST)
            || (refEntryType == vtkNewFoamToken::SCALARLIST
            && entry.FirstValue().GetType() == vtkNewFoamToken::VECTORLIST
            && nRefComponents == nComponents))
          {
          for (int tupleI = 0; tupleI < nTuples; tupleI++)
            {
            float *tuple = data->GetPointer(nComponents * tupleI);
            for (int componentI = 0; componentI < nComponents; componentI++)
              {
              tuple[componentI] += refTuple[componentI];
              }
            }
          }
        else
          {
          vtkErrorMacro(<< "Wrong referenceLevel type");
          return NULL;
          }
        }

#if vtksys_DATE_STAMP_FULL >= 20080620
      // swap the components of symmTensor to match the component
      // names in paraview
      const int nComponents = data->GetNumberOfComponents();
      if (nComponents == 6)
        {
        for (int tupleI = 0; tupleI < nTuples; tupleI++)
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
    else if (entry.FirstValue().GetType() == vtkNewFoamToken::EMPTYLIST && nElements <= 0)
      {
      data = vtkFloatArray::New();
      // set the number of components as appropriate if the list is empty
      if (fieldType == "ScalarField" || fieldType == "SphericalTensorField")
        {
        data->SetNumberOfComponents(1);
        }
      else if (fieldType == "VectorField")
        {
        data->SetNumberOfComponents(3);
        }
      else if (fieldType == "SymmTensorField")
        {
        data->SetNumberOfComponents(6);
        }
      else if (fieldType == "TensorField")
        {
        data->SetNumberOfComponents(9);
        }
      }
    else
      {
      vtkErrorMacro(<< ioPtr->GetFileName().c_str() << " is not a valid "
          << ioPtr->GetClassName().c_str());
      return NULL;
      }
    }
  else // undefined
    {
    vtkErrorMacro(<< "Either uniform or nonuniform keyword is missing in entry "
        << entry.GetKeyword().c_str() << " in " << ioPtr->GetFileName().c_str());
    return NULL;
    }
  return data;
}

//-----------------------------------------------------------------------------
// convert OpenFOAM's dimension array representation to string
void vtkNewOpenFOAMReaderPrivate::ConstructDimensions(vtkStdString *dimString,
    vtkNewFoamDict *dictPtr)
{
  if (!this->Parent->GetAddDimensionsToArrayNames())
    {
    return;
    }
  vtkNewFoamEntry *dimEntry = dictPtr->Lookup("dimensions");
  if (dimEntry != NULL && dimEntry->FirstValue().GetType() == vtkNewFoamToken::LABELLIST)
    {
    vtkIntArray &dims = dimEntry->LabelList();
    if (dims.GetNumberOfTuples() == 7)
      {
      int dimSet[7];
      for (int dimI = 0; dimI < 7; dimI++)
        {
        dimSet[dimI] = dims.GetValue(dimI);
        }
      static const char *units[7] =
      { "kg", "m", "s", "K", "mol", "A", "cd" };
      vtksys_ios::ostringstream posDim, negDim;
      int posSpc = 0, negSpc = 0;
      if (dimSet[0] == 1 && dimSet[1] == -1 && dimSet[2] == -2)
        {
        posDim << "Pa";
        dimSet[0] = dimSet[1] = dimSet[2] = 0;
        posSpc = 1;
        }
      for (int dimI = 0; dimI < 7; dimI++)
        {
        const int dimDim = dimSet[dimI];
        if (dimDim > 0)
          {
          if (posSpc)
            {
            posDim << " ";
            }
          posDim << units[dimI];
          if (dimDim > 1)
            {
            posDim << dimDim;
            }
          posSpc++;
          }
        else if (dimDim < 0)
          {
          if (negSpc)
            {
            negDim << " ";
            }
          negDim << units[dimI];
          if (dimDim < -1)
            {
            negDim << -dimDim;
            }
          negSpc++;
          }
        }
      *dimString += " [" + posDim.str();
      if (negSpc > 0)
        {
        if (posSpc == 0)
          {
          *dimString += "1";
          }
        if (negSpc > 1)
          {
          *dimString += "/(" + negDim.str() + ")";
          }
        else
          {
          *dimString += "/" + negDim.str();
          }
        }
      else if (posSpc == 0)
        {
        *dimString += "-";
        }
      *dimString += "]";
      }
    }
}

//-----------------------------------------------------------------------------
void vtkNewOpenFOAMReaderPrivate::GetVolFieldAtTimeStep(
    vtkUnstructuredGrid *internalMesh, vtkMultiBlockDataSet *boundaryMesh,
    const vtkStdString &varName)
{
  vtkNewFoamIOobject io(this->CasePath,
      this->Parent->GetIsSinglePrecisionBinary() != 0);
  vtkNewFoamDict dict;
  if (!this->ReadFieldFile(&io, &dict, varName,
      this->Parent->CellDataArraySelection))
    {
    return;
    }

  if (io.GetClassName().substr(0, 3) != "vol")
    {
    vtkErrorMacro(<< io.GetFileName().c_str() << " is not a volField");
    return;
    }

  vtkNewFoamEntry *iEntry = dict.Lookup("internalField");
  if (iEntry == NULL)
    {
    vtkErrorMacro(<<"internalField not found in " << io.GetFileName().c_str());
    return;
    }

  if (iEntry->FirstValue().GetType() == vtkNewFoamToken::EMPTYLIST)
    {
    // if there's no cell there shouldn't be any boundary faces either
    if (this->NumCells > 0)
      {
      vtkErrorMacro(<<"internalField of " << io.GetFileName().c_str()
          << " is empty");
      }
    return;
    }

  vtkNewFoamEntry *rEntry = dict.Lookup("referenceLevel");
  vtkStdString fieldType = io.GetClassName().substr(3, vtkStdString::npos);
  vtkFloatArray *iData =
      this->FillField(iEntry, rEntry, this->NumCells, &io, fieldType, false);
  if (iData == NULL)
    {
    return;
    }

  vtkStdString dimString;
  this->ConstructDimensions(&dimString, &dict);

  vtkFloatArray *ctpData = NULL;

  if (iData->GetSize() > 0)
    {
    // Add field only if internal Mesh exists (skip if not selected).
    // Note we still need to read internalField even if internal mesh is
    // not selected, since boundaries without value entries may refer to
    // the internalField.
    if (internalMesh != NULL)
      {
      if (this->Parent->GetDecomposePolyhedra())
        {
        // add values for decomposed cells
        this->ExtendArray<vtkFloatArray, float>(iData, this->NumCells
            + this->NumTotalAdditionalCells);
        const int nTuples = this->AdditionalCellIds->GetNumberOfTuples();
        int additionalCellI = this->NumCells;
        for (int tupleI = 0; tupleI < nTuples; tupleI++)
          {
          const int nCells = this->NumAdditionalCells->GetValue(tupleI);
          const vtkIdType cellId = this->AdditionalCellIds->GetValue(tupleI);
          for (int cellI = 0; cellI < nCells; cellI++)
            {
            iData->InsertTuple(additionalCellI++, cellId, iData);
            }
          }
        }

      // set data to internal mesh
      this->AddArrayToFieldData(internalMesh->GetCellData(), iData,
          io.GetObjectName() + dimString);

      if (this->Parent->GetCreateCellToPoint())
        {
        // Create cell-to-point interpolated data
        ctpData = vtkFloatArray::New();
        ctpData->SetNumberOfComponents(iData->GetNumberOfComponents());
        ctpData->SetNumberOfTuples(internalMesh->GetPoints()->GetNumberOfPoints());
        if (this->InternalPoints != NULL)
          {
          this->InterpolateCellToPoint(ctpData, iData, internalMesh,
              this->InternalPoints, this->InternalPoints->GetNumberOfTuples());
          }

        if (this->Parent->GetDecomposePolyhedra())
          {
          // assign cell values to additional points
          const int nPoints = this->AdditionalCellIds->GetNumberOfTuples();
          for (int pointI = 0; pointI < nPoints; pointI++)
            {
            ctpData->SetTuple(this->NumPoints + pointI,
                this->AdditionalCellIds->GetValue(pointI), iData);
            }
          }
        }
      }
    }
  else
    {
    // determine as there's no cells
    iData->Delete();
    return;
    }

  if (boundaryMesh == NULL)
    {
    iData->Delete();
    if (ctpData != NULL)
      {
      ctpData->Delete();
      }
    return;
    }

  vtkFloatArray *acData = NULL;

  if (this->Parent->GetCreateCellToPoint())
    {
    if (this->AllBoundaries == NULL)
      {
      vtkErrorMacro(<<"boundary mesh for cell to point filtering not found");
      iData->Delete();
      if (ctpData != NULL)
        {
        ctpData->Delete();
        }
      return;
      }
    acData = vtkFloatArray::New();
    acData->SetNumberOfComponents(iData->GetNumberOfComponents());
    acData->SetNumberOfTuples(this->AllBoundaries->GetNumberOfCells());
    }

  // set boundary values
  const vtkNewFoamEntry *bEntry = dict.Lookup("boundaryField");
  if (bEntry == NULL)
    {
    vtkErrorMacro(<< "boundaryField not found in object " << varName.c_str()
        << " at time = " << this->TimeNames->GetValue(this->TimeStep).c_str());
    iData->Delete();
    if (acData != NULL)
      {
      acData->Delete();
      }
    if (ctpData != NULL)
      {
      ctpData->Delete();
      }
    return;
    }

  vtkstd::vector<vtkFloatArray *> vDataVector;
  for (int boundaryI = 0, activeBoundaryI = 0; boundaryI
    < static_cast<int>(this->BoundaryDict.size()); boundaryI++)
    {
    const vtkNewFoamBoundaryEntry &beI = this->BoundaryDict[boundaryI];
    const vtkStdString &boundaryNameI = beI.BoundaryName;

    const vtkNewFoamEntry *bEntryI = bEntry->Dictionary().Lookup(boundaryNameI);
    if (bEntryI == NULL)
      {
      vtkErrorMacro(<< "boundaryField " << boundaryNameI.c_str()
          << " not found in object " << varName.c_str() << " at time = "
          << this->TimeNames->GetValue(this->TimeStep).c_str());
      iData->Delete();
      if (acData != NULL)
        {
        acData->Delete();
        }
      if (ctpData != NULL)
        {
        ctpData->Delete();
        }
      return;
      }

    if (bEntryI->FirstValue().GetType() != vtkNewFoamToken::DICTIONARY)
      {
      vtkErrorMacro(<< "Type of boundaryField " << boundaryNameI.c_str()
          << " is not a subdictionary in object " << varName.c_str()
          << " at time = " << this->TimeNames->GetValue(this->TimeStep).c_str());
      iData->Delete();
      if (acData != NULL)
        {
        acData->Delete();
        }
      if (ctpData != NULL)
        {
        ctpData->Delete();
        }
      return;
      }

    const int nFaces = beI.NFaces;

    vtkFloatArray* vData = NULL;
    if (beI.BoundaryType != vtkNewFoamBoundaryEntry::GEOMETRICAL)
      {
      vtkNewFoamEntry *vEntry = bEntryI->Dictionary().Lookup("value");
      if (vEntry != NULL) // the boundary has a value entry
        {
        vData = this->FillField(vEntry, rEntry, nFaces, &io, fieldType, false);
        if (vData == NULL)
          {
          iData->Delete();
          if (acData != NULL)
            {
            acData->Delete();
            }
          if (ctpData != NULL)
            {
            ctpData->Delete();
            }
          return;
          }
        }
      else
        {
        // uniformFixedValue B.C.
        const vtkNewFoamEntry *ufvEntry = bEntryI->Dictionary().Lookup("type");
        if (ufvEntry != NULL)
          {
          if (ufvEntry->ToWord() == "uniformFixedValue")
            {
            // the boundary is of uniformFixedValue type
            vtkNewFoamEntry *uvEntry = bEntryI->Dictionary().Lookup("uniformValue");
            if (uvEntry != NULL) // and has a uniformValue entry
              {
              vData = this->FillField(uvEntry, rEntry, nFaces, &io, fieldType, true);
              if (vData == NULL)
                {
                iData->Delete();
                if (acData != NULL)
                  {
                  acData->Delete();
                  }
                if (ctpData != NULL)
                  {
                  ctpData->Delete();
                  }
                return;
                }
              }
            }
          }
        }
      }
    const int boundaryStartFace = beI.StartFace
        - this->BoundaryDict[0].StartFace;

    if (vData == NULL) // doesn't have a value nor uniformValue entry
      {
      // use patch-internal values as boundary values
      vData = vtkFloatArray::New();
      vData->SetNumberOfComponents(iData->GetNumberOfComponents());
      vData->SetNumberOfTuples(nFaces);
      for (int j = 0; j < nFaces; j++)
        {
        const int cellId = this->FaceOwner->GetValue(boundaryStartFace + j);
        vData->SetTuple(j, cellId, iData);
        }
      }

    if (this->Parent->GetCreateCellToPoint())
      {
      const int startFace = beI.AllBoundariesStartFace;
      // if reading a processor sub-case of a decomposed case as is,
      // use the patch values of the processor patch as is
      if (beI.BoundaryType == vtkNewFoamBoundaryEntry::PHYSICAL
          || (this->ProcessorName == "" && beI.BoundaryType
              == vtkNewFoamBoundaryEntry::PROCESSOR))
        {
        // set the same value to AllBoundaries
        for (int faceI = 0; faceI < nFaces; faceI++)
          {
          acData->SetTuple(faceI + startFace, faceI, vData);
          }
        }
      // implies && this->ProcessorName != ""
      else if (beI.BoundaryType == vtkNewFoamBoundaryEntry::PROCESSOR)
        {
        // average patch internal value and patch value assuming the
        // patch value to be the patchInternalField of the neighbor
        // decomposed mesh. Using double precision to avoid degrade in
        // accuracy.
        const int nComponents = vData->GetNumberOfComponents();
        for (int faceI = 0; faceI < nFaces; faceI++)
          {
          const float *vTuple = vData->GetPointer(nComponents * faceI);
          const float *iTuple = iData->GetPointer(nComponents
              * this->FaceOwner->GetValue(boundaryStartFace + faceI));
          float *acTuple =
              acData->GetPointer(nComponents * (startFace + faceI));
          for (int componentI = 0; componentI < nComponents; componentI++)
            {
            acTuple[componentI] = (static_cast<double>(vTuple[componentI])
                + static_cast<double>(iTuple[componentI])) * 0.5;
            }
          }
        }
      }

    if (beI.IsActive)
      {
      vtkPolyData *bm =
          vtkPolyData::SafeDownCast(boundaryMesh->GetBlock(activeBoundaryI));
      this->AddArrayToFieldData(bm->GetCellData(), vData, io.GetObjectName()
          + dimString);

      if (this->Parent->GetCreateCellToPoint())
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
        this->AddArrayToFieldData(bm->GetPointData(), pData, io.GetObjectName()
            + dimString);
        pData->Delete();
        }

      activeBoundaryI++;
      }
    vData->Delete();
    }
  iData->Delete();

  if (this->Parent->GetCreateCellToPoint())
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

    if (ctpData != NULL)
      {
      // set cell-to-pint data for internal mesh
      for (int pointI = 0; pointI < nPoints; pointI++)
        {
        ctpData->SetTuple(this->AllBoundariesPointMap->GetValue(pointI),
            pointI, bpData);
        }
      this->AddArrayToFieldData(internalMesh->GetPointData(), ctpData,
          io.GetObjectName() + dimString);
      ctpData->Delete();
      }

    bpData->Delete();
    }
}

//-----------------------------------------------------------------------------
// read point field at a timestep
void vtkNewOpenFOAMReaderPrivate::GetPointFieldAtTimeStep(
    vtkUnstructuredGrid *internalMesh, vtkMultiBlockDataSet *boundaryMesh,
    const vtkStdString &varName)
{
  vtkNewFoamIOobject io(this->CasePath,
      this->Parent->GetIsSinglePrecisionBinary() != 0);
  vtkNewFoamDict dict;
  if (!this->ReadFieldFile(&io, &dict, varName,
      this->Parent->PointDataArraySelection))
    {
    return;
    }

  if (io.GetClassName().substr(0, 5) != "point")
    {
    vtkErrorMacro(<< io.GetFileName().c_str() << " is not a pointField");
    return;
    }

  vtkNewFoamEntry *iEntry = dict.Lookup("internalField");
  if (iEntry == NULL)
    {
    vtkErrorMacro(<<"internalField not found in " << io.GetFileName().c_str());
    return;
    }

  if (iEntry->FirstValue().GetType() == vtkNewFoamToken::EMPTYLIST)
    {
    // if there's no cell there shouldn't be any boundary faces either
    if (this->NumPoints > 0)
      {
      vtkErrorMacro(<<"internalField of " << io.GetFileName().c_str()
          << " is empty");
      }
    return;
    }

  vtkNewFoamEntry *rEntry = dict.Lookup("referenceLevel");
  vtkStdString fieldType = io.GetClassName().substr(5, vtkStdString::npos);
  vtkFloatArray *iData = this->FillField(iEntry, rEntry, this->NumPoints, &io,
      fieldType, false);
  if (iData == NULL)
    {
    return;
    }

  vtkStdString dimString;
  this->ConstructDimensions(&dimString, &dict);

  // AdditionalCellPoints is NULL if creation of InternalMesh had been skipped
  if (this->AdditionalCellPoints != NULL)
    {
    // point-to-cell interpolation to additional cell centroidal points
    // for decomposed cells
    const int nAdditionalPoints = static_cast<int>(this->AdditionalCellPoints->size());
    const int nComponents = iData->GetNumberOfComponents();
    this->ExtendArray<vtkFloatArray, float>(iData, this->NumPoints
        + nAdditionalPoints);
    for (int i = 0; i < nAdditionalPoints; i++)
      {
      vtkIntArray *acp = this->AdditionalCellPoints->operator[](i);
      int nPoints = acp->GetDataSize();
      double interpolatedValue[9];
      for (int k = 0; k < nComponents; k++)
        {
        interpolatedValue[k] = 0.0;
        }
      for (int j = 0; j < nPoints; j++)
        {
        const float *tuple = iData->GetPointer(nComponents * acp->GetValue(j));
        for (int k = 0; k < nComponents; k++)
          {
          interpolatedValue[k] += tuple[k];
          }
        }
      const double weight = 1.0 / static_cast<double>(nPoints);
      for (int k = 0; k < nComponents; k++)
        {
        interpolatedValue[k] *= weight;
        }
      // will automatically be converted to float
      iData->InsertTuple(this->NumPoints + i, interpolatedValue);
      }
    }

  if (iData->GetSize() > 0)
    {
    // Add field only if internal Mesh exists (skip if not selected).
    // Note we still need to read internalField even if internal mesh is
    // not selected, since boundaries without value entries may refer to
    // the internalField.
    if (internalMesh != NULL)
      {
      // set data to internal mesh
      this->AddArrayToFieldData(internalMesh->GetPointData(), iData,
          io.GetObjectName() + dimString);
      }
    }
  else
    {
    // determine as there's no points
    iData->Delete();
    return;
    }

  if (boundaryMesh == NULL)
    {
    iData->Delete();
    return;
    }

  // use patch-internal values as boundary values
  for (vtkIdType boundaryI = 0, activeBoundaryI = 0; boundaryI
    < static_cast<vtkIdType>(this->BoundaryDict.size()); boundaryI++)
    {
    if (this->BoundaryDict[boundaryI].IsActive)
      {
      vtkFloatArray *vData = vtkFloatArray::New();
      vtkIntArray& bpMap = *this->BoundaryPointMap->operator[](activeBoundaryI);
      const int nPoints = bpMap.GetNumberOfTuples();
      vData->SetNumberOfComponents(iData->GetNumberOfComponents());
      vData->SetNumberOfTuples(nPoints);
      for (int j = 0; j < nPoints; j++)
        {
        vData->SetTuple(j, bpMap.GetValue(j), iData);
        }
      this->AddArrayToFieldData(vtkPolyData::SafeDownCast(
          boundaryMesh->GetBlock(activeBoundaryI))->GetPointData(), vData, io.GetObjectName()
          + dimString);
      vData->Delete();
      activeBoundaryI++;
      }
    }
  iData->Delete();
}

//-----------------------------------------------------------------------------
vtkMultiBlockDataSet* vtkNewOpenFOAMReaderPrivate::MakeLagrangianMesh()
{
  vtkMultiBlockDataSet *lagrangianMesh = vtkMultiBlockDataSet::New();

  for (int cloudI = 0; cloudI
      < this->Parent->LagrangianPaths->GetNumberOfTuples(); cloudI++)
    {
    const vtkStdString& pathI = this->Parent->LagrangianPaths->GetValue(cloudI);

    // still can't distinguish on patch selection panel, but can
    // distinguish the "lagrangian" reserved path component and a mesh
    // region with the same name
    vtkStdString subCloudName;
    if (pathI[0] == '/')
      {
      subCloudName = pathI.substr(1, vtkStdString::npos);
      }
    else
      {
      subCloudName = pathI;
      }
    if (this->RegionName != pathI.substr(0, pathI.find('/'))
        || !this->Parent->GetPatchArrayStatus(subCloudName.c_str()))
      {
      continue;
      }

    const vtkStdString cloudPath(this->CurrentTimePath() + "/" + subCloudName
        + "/");
    const vtkStdString positionsPath(cloudPath + "positions");

    // create an empty mesh to keep node/leaf structure of the
    // multi-block consistent even if mesh doesn't exist
    vtkPolyData *meshI = vtkPolyData::New();
    const int blockI = lagrangianMesh->GetNumberOfBlocks();
    lagrangianMesh->SetBlock(blockI, meshI);
    // extract the cloud name
    this->SetBlockName(lagrangianMesh, blockI, pathI.substr(pathI.rfind('/') + 1).c_str());

    vtkNewFoamIOobject io(this->CasePath,
        this->Parent->GetIsSinglePrecisionBinary() != 0);
    if (!(io.Open(positionsPath) || io.Open(positionsPath + ".gz")))
      {
      meshI->Delete();
      continue;
      }

    // tell the IO object if the file is in OF 1.3 binary
    // lagrangian/positions format
    io.SetIs13Positions(this->Parent->GetPositionsIsIn13Format() != 0);

    vtkNewFoamEntryValue dict(NULL);
    try
      {
      dict.ReadNonuniformList<vtkNewFoamToken::VECTORLIST,
      vtkNewFoamEntryValue::vectorListTraits<vtkFloatArray, float, 3, true> >(
          io);
      }
    catch(vtkNewFoamError& e)
      {
      vtkErrorMacro(<<"Error reading line " << io.GetLineNumber()
          << " of " << io.GetFileName().c_str() << ": " << e.c_str());
      meshI->Delete();
      continue;
      }
    io.Close();

    vtkFloatArray *pointArray = static_cast<vtkFloatArray *>(dict.Ptr());
    const int nParticles = pointArray->GetNumberOfTuples();

    // instantiate the points class
    vtkPoints *points = vtkPoints::New();
    points->SetData(pointArray);
    pointArray->Delete();

    // create lagrangian mesh
    meshI->Allocate(nParticles);
    for (vtkIdType i = 0; i < nParticles; i++)
      {
      meshI->InsertNextCell(VTK_VERTEX, 1, &i);
      }
    meshI->SetPoints(points);
    points->Delete();

    // read lagrangian fields
    for (int fieldI = 0; fieldI
        < this->LagrangianFieldFiles->GetNumberOfValues(); fieldI++)
      {
      const vtkStdString varPath(cloudPath
          + this->LagrangianFieldFiles->GetValue(fieldI));

      vtkNewFoamIOobject io2(this->CasePath,
          this->Parent->GetIsSinglePrecisionBinary() != 0);
      if (!io2.Open(varPath))
        {
        // if the field file doesn't exist we simply return without
        // issuing an error as a simple way of supporting multi-region
        // lagrangians
        continue;
        }

      // if the variable is disabled on selection panel then skip it
      const vtkStdString selectionName(io2.GetObjectName());
      if (this->Parent->LagrangianDataArraySelection->ArrayExists(selectionName.c_str())
          && !this->Parent->GetLagrangianArrayStatus(selectionName.c_str()))
        {
        continue;
        }

      // read the field file into dictionary
      vtkNewFoamEntryValue dict2(NULL);
      if (!dict2.ReadField(io2))
        {
        vtkErrorMacro(<<"Error reading line " << io2.GetLineNumber()
            << " of " << io2.GetFileName().c_str() << ": "
            << io2.GetError().c_str());
        continue;
        }

      // set lagrangian values
      if (dict2.GetType() != vtkNewFoamToken::SCALARLIST && dict2.GetType()
          != vtkNewFoamToken::VECTORLIST && dict2.GetType()
          != vtkNewFoamToken::LABELLIST)
        {
        vtkErrorMacro(<< io2.GetFileName().c_str()
            << ": Unsupported lagrangian field type "
            << io2.GetClassName().c_str());
        continue;
        }

      vtkDataArray* lData = static_cast<vtkDataArray *>(dict2.Ptr());

      // GetNumberOfTuples() works for both scalar and vector
      const int nParticles2 = lData->GetNumberOfTuples();
      if (nParticles2 != meshI->GetNumberOfCells())
        {
        vtkErrorMacro(<< io2.GetFileName().c_str()
            <<": Sizes of lagrangian mesh and field don't match: mesh = "
            << meshI->GetNumberOfCells() << ", field = " << nParticles2);
        lData->Delete();
        continue;
        }

      this->AddArrayToFieldData(meshI->GetCellData(), lData, selectionName);
      if (this->Parent->GetCreateCellToPoint())
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
vtkNewFoamDict* vtkNewOpenFOAMReaderPrivate::GatherBlocks(const char* typeIn,
    const int timeStep)
{
  if (this->PolyMeshFacesDir->GetValue(timeStep) == "")
    {
    return NULL;
    }

  vtkStdString type(typeIn);
  vtkStdString blockPath =
      this->TimeRegionMeshPath(this->PolyMeshFacesDir, timeStep) + type;

  vtkNewFoamIOobject io(this->CasePath,
      this->Parent->GetIsSinglePrecisionBinary() != 0);
  if (!(io.Open(blockPath) || io.Open(blockPath + ".gz")))
    {
    return NULL;
    }

  vtkNewFoamDict* dictPtr = new vtkNewFoamDict;
  vtkNewFoamDict& dict = *dictPtr;
  if (!dict.Read(io))
    {
    vtkErrorMacro(<<"Error reading line " << io.GetLineNumber()
        << " of " << io.GetFileName().c_str() << ": " << io.GetError().c_str());
    delete dictPtr;
    return NULL;
    }
  if (dict.GetType() != vtkNewFoamToken::DICTIONARY)
    {
    vtkErrorMacro(<<"The file type of " << io.GetFileName().c_str()
        << " is not a dictionary");
    delete dictPtr;
    return NULL;
    }
  return dictPtr;
}

//-----------------------------------------------------------------------------
// returns a requested point zone mesh
bool vtkNewOpenFOAMReaderPrivate::GetPointZoneMesh(
    vtkMultiBlockDataSet *pointZoneMesh, vtkPoints *points)
{
  vtkNewFoamDict *pointZoneDictPtr
      = this->GatherBlocks("pointZones", this->TimeStep);

  if (pointZoneDictPtr == NULL)
    {
    // not an error
    return true;
    }

  vtkNewFoamDict &pointZoneDict = *pointZoneDictPtr;
  int nPointZones = static_cast<int>(pointZoneDict.size());

  for (int i = 0; i < nPointZones; i++)
    {
    // look up point labels
    vtkNewFoamDict &dict = pointZoneDict[i]->Dictionary();
    vtkNewFoamEntry *pointLabelsEntry = dict.Lookup("pointLabels");
    if (pointLabelsEntry == NULL)
      {
      delete pointZoneDictPtr;
      vtkErrorMacro(<<"pointLabels not found in pointZones");
      return false;
      }

    // allocate an empty mesh if the list is empty
    if (pointLabelsEntry->FirstValue().GetType() == vtkNewFoamToken::EMPTYLIST)
      {
      vtkPolyData *pzm = vtkPolyData::New();
      pointZoneMesh->SetBlock(i, pzm);
      pzm->Delete();
      // set name
      this->SetBlockName(pointZoneMesh, i, pointZoneDict[i]->GetKeyword().c_str());
      continue;
      }

    if (pointLabelsEntry->FirstValue().GetType() != vtkNewFoamToken::LABELLIST)
      {
      delete pointZoneDictPtr;
      vtkErrorMacro(<<"pointLabels not of type labelList: type = "
          << pointLabelsEntry->FirstValue().GetType());
      return false;
      }

    vtkIntArray &labels = pointLabelsEntry->LabelList();

    int nPoints = labels.GetNumberOfTuples();
    if (nPoints > this->NumPoints)
      {
      vtkErrorMacro(<<"The length of pointLabels " << nPoints
          << " for pointZone " << pointZoneDict[i]->GetKeyword().c_str()
          << " exceeds the number of points " << this->NumPoints);
      delete pointZoneDictPtr;
      return false;
      }

    // allocate new grid: we do not use resize() beforehand since it
    // could lead to undefined pointer if we return by error
    vtkPolyData *pzm = vtkPolyData::New();

    // set pointZone size
    pzm->Allocate(nPoints);

    // insert points
    for (int j = 0; j < nPoints; j++)
      {
      vtkIdType pointLabel = labels.GetValue(j); // must be vtkIdType
      if (pointLabel >= this->NumPoints)
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
    this->SetBlockName(pointZoneMesh, i, pointZoneDict[i]->GetKeyword().c_str());
    }

  delete pointZoneDictPtr;

  return true;
}

//-----------------------------------------------------------------------------
// returns a requested face zone mesh
bool vtkNewOpenFOAMReaderPrivate::GetFaceZoneMesh(
    vtkMultiBlockDataSet *faceZoneMesh,
    const vtkNewFoamIntVectorVector *facesPoints, vtkPoints *points)
{
  vtkNewFoamDict *faceZoneDictPtr
      = this->GatherBlocks("faceZones", this->TimeStep);

  if (faceZoneDictPtr == NULL)
    {
    // not an error
    return true;
    }

  vtkNewFoamDict &faceZoneDict = *faceZoneDictPtr;
  int nFaceZones = static_cast<int>(faceZoneDict.size());

  for (int i = 0; i < nFaceZones; i++)
    {
    // look up face labels
    vtkNewFoamDict &dict = faceZoneDict[i]->Dictionary();
    vtkNewFoamEntry *faceLabelsEntry = dict.Lookup("faceLabels");
    if (faceLabelsEntry == NULL)
      {
      delete faceZoneDictPtr;
      vtkErrorMacro(<<"faceLabels not found in faceZones");
      return false;
      }

    // allocate an empty mesh if the list is empty
    if (faceLabelsEntry->FirstValue().GetType() == vtkNewFoamToken::EMPTYLIST)
      {
      vtkPolyData *fzm = vtkPolyData::New();
      faceZoneMesh->SetBlock(i, fzm);
      fzm->Delete();
      // set name
      this->SetBlockName(faceZoneMesh, i, faceZoneDict[i]->GetKeyword().c_str());
      continue;
      }

    if (faceLabelsEntry->FirstValue().GetType() != vtkNewFoamToken::LABELLIST)
      {
      delete faceZoneDictPtr;
      vtkErrorMacro(<<"faceLabels not of type labelList");
      return false;
      }

    vtkIntArray &labels = faceLabelsEntry->LabelList();

    int nFaces = labels.GetNumberOfTuples();
    if (nFaces > this->FaceOwner->GetNumberOfTuples())
      {
      vtkErrorMacro(<<"The length of faceLabels " << nFaces
          << " for faceZone " << faceZoneDict[i]->GetKeyword().c_str()
          << " exceeds the number of faces "
          << this->FaceOwner->GetNumberOfTuples());
      delete faceZoneDictPtr;
      return false;
      }

    // allocate new grid: we do not use resize() beforehand since it
    // could lead to undefined pointer if we return by error
    vtkPolyData *fzm = vtkPolyData::New();

    // set faceZone size
    fzm->Allocate(nFaces);

    // aloocate array for converting int vector to vtkIdType vector:
    // workaround for 64bit machines
    int maxNFacePoints = 0;
    for (int j = 0; j < nFaces; j++)
      {
      const int nFacePoints = facesPoints->GetSize(labels.GetValue(j));
      if (nFacePoints > maxNFacePoints)
        {
        maxNFacePoints = nFacePoints;
        }
      }
    vtkIdList *facePointsVtkId = vtkIdList::New();
    facePointsVtkId->SetNumberOfIds(maxNFacePoints);

    // insert faces
    if (!this->InsertFacesToGrid(fzm, facesPoints, 0, nFaces, NULL,
        facePointsVtkId, &labels, false))
      {
      facePointsVtkId->Delete();
      fzm->Delete();
      return false;
      }
    facePointsVtkId->Delete();
    fzm->SetPoints(points);
    faceZoneMesh->SetBlock(i, fzm);
    fzm->Delete();
    // set name
    this->SetBlockName(faceZoneMesh, i, faceZoneDict[i]->GetKeyword().c_str());
    }

  delete faceZoneDictPtr;

  return true;
}

//-----------------------------------------------------------------------------
// returns a requested cell zone mesh
bool vtkNewOpenFOAMReaderPrivate::GetCellZoneMesh(
    vtkMultiBlockDataSet *cellZoneMesh,
    const vtkNewFoamIntVectorVector *cellsFaces,
    const vtkNewFoamIntVectorVector *facesPoints, vtkPoints *points)
{
  vtkNewFoamDict *cellZoneDictPtr
      = this->GatherBlocks("cellZones", this->TimeStep);

  if (cellZoneDictPtr == NULL)
    {
    // not an error
    return true;
    }

  vtkNewFoamDict &cellZoneDict = *cellZoneDictPtr;
  int nCellZones = static_cast<int>(cellZoneDict.size());

  for (int i = 0; i < nCellZones; i++)
    {
    // look up cell labels
    vtkNewFoamDict &dict = cellZoneDict[i]->Dictionary();
    vtkNewFoamEntry *cellLabelsEntry = dict.Lookup("cellLabels");
    if (cellLabelsEntry == NULL)
      {
      delete cellZoneDictPtr;
      vtkErrorMacro(<<"cellLabels not found in cellZones");
      return false;
      }

    // allocate an empty mesh if the list is empty
    if (cellLabelsEntry->FirstValue().GetType() == vtkNewFoamToken::EMPTYLIST)
      {
      vtkUnstructuredGrid *czm = vtkUnstructuredGrid::New();
      cellZoneMesh->SetBlock(i, czm);
      // set name
      this->SetBlockName(cellZoneMesh, i, cellZoneDict[i]->GetKeyword().c_str());
      continue;
      }

    if (cellLabelsEntry->FirstValue().GetType() != vtkNewFoamToken::LABELLIST)
      {
      delete cellZoneDictPtr;
      vtkErrorMacro(<<"cellLabels not of type labelList");
      return false;
      }

    vtkIntArray &labels = cellLabelsEntry->LabelList();

    int nCells = labels.GetNumberOfTuples();
    if (nCells > this->NumCells)
      {
      vtkErrorMacro(<<"The length of cellLabels " << nCells
          << " for cellZone " << cellZoneDict[i]->GetKeyword().c_str()
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
    this->SetBlockName(cellZoneMesh, i, cellZoneDict[i]->GetKeyword().c_str());
    }

  delete cellZoneDictPtr;
  return true;
}

//-----------------------------------------------------------------------------
void vtkNewOpenFOAMReaderPrivate::AddArrayToFieldData(
    vtkDataSetAttributes *fieldData, vtkDataArray *array,
    const vtkStdString &arrayName)
{
  // exclude dimensional unit string if any
  const vtkStdString arrayNameString(arrayName.substr(0, arrayName.find(' ')));
  array->SetName(arrayName.c_str());

  if (array->GetNumberOfComponents() == 1 && arrayNameString == "p")
    {
    fieldData->SetScalars(array);
    }
  else if (array->GetNumberOfComponents() == 3 && arrayNameString == "U")
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
int vtkNewOpenFOAMReaderPrivate::RequestData(vtkMultiBlockDataSet *output,
    bool recreateInternalMesh, bool recreateBoundaryMesh, bool updateVariables)
{
  recreateInternalMesh |= this->TimeStepOld == -1
      // the following two quite likely indicate reading mesh failed on
      // the previous RequestData() call, hence trying again
      || (this->InternalMeshSelectionStatus && this->InternalMesh == NULL)
      || this->BoundaryMesh == NULL
      || this->InternalMeshSelectionStatus
          != this->InternalMeshSelectionStatusOld
      || this->PolyMeshFacesDir->GetValue(this->TimeStep)
          != this->PolyMeshFacesDir->GetValue(this->TimeStepOld)
      || this->FaceOwner == NULL;
  recreateBoundaryMesh |= recreateInternalMesh;
  updateVariables |= recreateBoundaryMesh || this->TimeStep
      != this->TimeStepOld;
  const bool pointsMoved = this->TimeStepOld == -1
      || this->PolyMeshPointsDir->GetValue(this->TimeStep)
          != this->PolyMeshPointsDir->GetValue(this->TimeStepOld);
  const bool moveInternalPoints = !recreateInternalMesh && pointsMoved;
  const bool moveBoundaryPoints = !recreateBoundaryMesh && pointsMoved;

  const bool createEulerians
      = this->PolyMeshFacesDir->GetValue(this->TimeStep) != "";

  // determine if we need to reconstruct meshes
  if (recreateInternalMesh)
    {
    this->ClearInternalMeshes();
    }
  if (recreateBoundaryMesh)
    {
    this->ClearBoundaryMeshes();
    }

  vtkNewFoamIntVectorVector *facePoints = NULL;
  vtkStdString meshDir;
  if (createEulerians && (recreateInternalMesh || recreateBoundaryMesh))
    {
    // create paths to polyMesh files
    meshDir = this->CurrentTimeRegionMeshPath(this->PolyMeshFacesDir);

    // create the faces vector
    facePoints = this->ReadFacesFile(meshDir);
    if (facePoints == NULL)
      {
      return 0;
      }
    this->Parent->UpdateProgress(0.2);
    }

  vtkNewFoamIntVectorVector *cellFaces = NULL;
  if (createEulerians && recreateInternalMesh)
    {
    // read owner/neighbor and create the FaceOwner and cellFaces vectors
    cellFaces = this->ReadOwnerNeighborFiles(meshDir, facePoints);
    if (cellFaces == NULL)
      {
      delete facePoints;
      return 0;
      }
    this->Parent->UpdateProgress(0.3);
    }

  vtkFloatArray *pointArray = NULL;
  if (createEulerians && (recreateInternalMesh || (recreateBoundaryMesh
      && !recreateInternalMesh && this->InternalMesh == NULL)
      || moveInternalPoints || moveBoundaryPoints))
    {
    // get the points
    pointArray = this->ReadPointsFile();
    if ((pointArray == NULL && recreateInternalMesh) || (facePoints != NULL
        && !this->CheckFacePoints(facePoints)))
      {
      delete cellFaces;
      delete facePoints;
      return 0;
      }
    this->Parent->UpdateProgress(0.4);
    }

  // make internal mesh
  // Create Internal Mesh only if required for display
  if (createEulerians && recreateInternalMesh)
    {
    if (this->Parent->GetPatchArrayStatus((this->RegionPrefix() + "internalMesh").c_str()))
      {
      this->InternalMesh = this->MakeInternalMesh(cellFaces, facePoints,
          pointArray);
      }
    // read and construct zones
    if (this->Parent->GetReadZones())
      {
      vtkPoints *points;
      if (this->InternalMesh != NULL)
        {
        points = this->InternalMesh->GetPoints();
        }
      else
        {
        points = vtkPoints::New();
        points->SetData(pointArray);
        }

      this->PointZoneMesh = vtkMultiBlockDataSet::New();
      if (!this->GetPointZoneMesh(this->PointZoneMesh, points))
        {
        this->PointZoneMesh->Delete();
        this->PointZoneMesh = NULL;
        delete cellFaces;
        delete facePoints;
        if (this->InternalMesh == NULL)
          {
          points->Delete();
          }
        pointArray->Delete();
        return 0;
        }
      if (this->PointZoneMesh->GetNumberOfBlocks() == 0)
        {
        this->PointZoneMesh->Delete();
        this->PointZoneMesh = NULL;
        }

      this->FaceZoneMesh = vtkMultiBlockDataSet::New();
      if (!this->GetFaceZoneMesh(this->FaceZoneMesh, facePoints, points))
        {
        this->FaceZoneMesh->Delete();
        this->FaceZoneMesh = NULL;
        if (this->PointZoneMesh != NULL)
          {
          this->PointZoneMesh->Delete();
          this->PointZoneMesh = NULL;
          }
        delete cellFaces;
        delete facePoints;
        if (this->InternalMesh == NULL)
          {
          points->Delete();
          }
        pointArray->Delete();
        return 0;
        }
      if (this->FaceZoneMesh->GetNumberOfBlocks() == 0)
        {
        this->FaceZoneMesh->Delete();
        this->FaceZoneMesh = NULL;
        }

      this->CellZoneMesh = vtkMultiBlockDataSet::New();
      if (!this->GetCellZoneMesh(this->CellZoneMesh, cellFaces, facePoints,
          points))
        {
        this->CellZoneMesh->Delete();
        this->CellZoneMesh = NULL;
        if (this->FaceZoneMesh != NULL)
          {
          this->FaceZoneMesh->Delete();
          this->FaceZoneMesh = NULL;
          }
        if (this->PointZoneMesh != NULL)
          {
          this->PointZoneMesh->Delete();
          this->PointZoneMesh = NULL;
          }
        delete cellFaces;
        delete facePoints;
        if (this->InternalMesh == NULL)
          {
          points->Delete();
          }
        pointArray->Delete();
        return 0;
        }
      if (this->CellZoneMesh->GetNumberOfBlocks() == 0)
        {
        this->CellZoneMesh->Delete();
        this->CellZoneMesh = NULL;
        }
      if (this->InternalMesh == NULL)
        {
        points->Delete();
        }
      }
    delete cellFaces;
    this->TruncateFaceOwner();
    }

  if (createEulerians && recreateBoundaryMesh)
    {
    vtkFloatArray *boundaryPointArray;
    if (pointArray != NULL)
      {
      boundaryPointArray = pointArray;
      }
    else
      {
      boundaryPointArray
          = static_cast<vtkFloatArray *>(this->InternalMesh->GetPoints()->GetData());
      }
    // create boundary mesh
    this->BoundaryMesh = this->MakeBoundaryMesh(facePoints, boundaryPointArray);
    if (this->BoundaryMesh == NULL)
      {
      delete facePoints;
      if (pointArray != NULL)
        {
        pointArray->Delete();
        }
      return 0;
      }
#if 0
    this->CalculateReciprocalDelta(facePoints);
#endif
    }

  delete facePoints;

  // if only point coordinates change refresh point vector
  if (createEulerians && moveInternalPoints)
    {
    // refresh the points in each mesh
    vtkPoints *points;
    // Check if Internal Mesh exists first....
    if (this->InternalMesh != NULL)
      {
      points = this->MoveInternalMesh(this->InternalMesh, pointArray);
      if (points == NULL)
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

    if (this->PointZoneMesh != NULL)
      {
      for (unsigned int i = 0; i < this->PointZoneMesh->GetNumberOfBlocks(); i++)
        {
        vtkPolyData::SafeDownCast(this->PointZoneMesh->GetBlock(i))
        ->SetPoints(points);
        }
      }
    if (this->FaceZoneMesh != NULL)
      {
      for (unsigned int i = 0; i < this->FaceZoneMesh->GetNumberOfBlocks(); i++)
        {
        vtkPolyData::SafeDownCast(this->FaceZoneMesh->GetBlock(i))
        ->SetPoints(points);
        }
      }
    if (this->CellZoneMesh != NULL)
      {
      for (unsigned int i = 0; i < this->CellZoneMesh->GetNumberOfBlocks(); i++)
        {
        vtkUnstructuredGrid::SafeDownCast(this->CellZoneMesh->GetBlock(i))
        ->SetPoints(points);
        }
      }
    points->Delete();
    }

  if (createEulerians && moveBoundaryPoints)
    {
    // Check if Boundary Mesh exists first....
    if (this->BoundaryMesh != NULL)
      {
      this->MoveBoundaryMesh(this->BoundaryMesh, pointArray);
      }
    }

  if (pointArray != NULL)
    {
    pointArray->Delete();
    }
  this->Parent->UpdateProgress(0.5);

  vtkMultiBlockDataSet *lagrangianMesh = NULL;
  if (updateVariables)
    {
    if (createEulerians)
      {
      if (!recreateInternalMesh && this->InternalMesh != NULL)
        {
        // clean up arrays of the previous timestep
        // Check if Internal Mesh Exists first...
        this->InternalMesh->GetCellData()->Initialize();
        this->InternalMesh->GetPointData()->Initialize();
        }
      // Check if Boundary Mesh Exists first...
      if (!recreateBoundaryMesh && this->BoundaryMesh != NULL)
        {
        for (unsigned int i = 0; i < this->BoundaryMesh->GetNumberOfBlocks(); i++)
          {
          vtkPolyData *bm =
              vtkPolyData::SafeDownCast(this->BoundaryMesh->GetBlock(i));
          bm->GetCellData()->Initialize();
          bm->GetPointData()->Initialize();
          }
        }
      // read field data variables into Internal/Boundary meshes
      for (int i = 0; i < (int)this->VolFieldFiles->GetNumberOfValues(); i++)
        {
        this->GetVolFieldAtTimeStep(this->InternalMesh, this->BoundaryMesh,
            this->VolFieldFiles->GetValue(i));
        this->Parent->UpdateProgress(0.5 + 0.25 * ((float)(i + 1)
            / ((float)this->VolFieldFiles->GetNumberOfValues() + 0.0001)));
        }
      for (int i = 0; i < (int)this->PointFieldFiles->GetNumberOfValues(); i++)
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
  if (this->InternalMesh != NULL)
    {
    output->SetBlock(0, this->InternalMesh);
    this->SetBlockName(output, 0, "internalMesh");
    }

  // set boundary meshes/data as output
  if (this->BoundaryMesh != NULL && this->BoundaryMesh->GetNumberOfBlocks() > 0)
    {
    const unsigned int groupTypeI = output->GetNumberOfBlocks();
    output->SetBlock(groupTypeI, this->BoundaryMesh);
    this->SetBlockName(output, groupTypeI, "Patches");
    }

  // set lagrangian mesh as output
  if (lagrangianMesh != NULL)
    {
    if (lagrangianMesh->GetNumberOfBlocks() > 0)
      {
      const unsigned int groupTypeI = output->GetNumberOfBlocks();
      output->SetBlock(groupTypeI, lagrangianMesh);
      this->SetBlockName(output, groupTypeI, "Lagrangian Particles");
      }
    lagrangianMesh->Delete();
    }

  if (this->Parent->GetReadZones())
    {
    vtkMultiBlockDataSet *zones = NULL;
    // set Zone Meshes as output
    if (this->PointZoneMesh != NULL)
      {
      zones = vtkMultiBlockDataSet::New();
      const unsigned int zoneTypeI = zones->GetNumberOfBlocks();
      zones->SetBlock(zoneTypeI, this->PointZoneMesh);
      this->SetBlockName(zones, zoneTypeI, "pointZones");
      }

    if (this->FaceZoneMesh != NULL)
      {
      if (zones == NULL)
        {
        zones = vtkMultiBlockDataSet::New();
        }
      const unsigned int zoneTypeI = zones->GetNumberOfBlocks();
      zones->SetBlock(zoneTypeI, this->FaceZoneMesh);
      this->SetBlockName(zones, zoneTypeI, "faceZones");
      }

    if (this->CellZoneMesh != NULL)
      {
      if (zones == NULL)
        {
        zones = vtkMultiBlockDataSet::New();
        }
      const unsigned int zoneTypeI = zones->GetNumberOfBlocks();
      zones->SetBlock(zoneTypeI, this->CellZoneMesh);
      this->SetBlockName(zones, zoneTypeI, "cellZones");
      }
    if (zones != NULL)
      {
      const unsigned int groupTypeI = output->GetNumberOfBlocks();
      output->SetBlock(groupTypeI, zones);
      this->SetBlockName(output, groupTypeI, "Zones");
      }
    }

  if (this->Parent->GetCacheMesh())
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
vtkNewOpenFOAMReader::vtkNewOpenFOAMReader()
{
  this->SetNumberOfInputPorts(0);

  this->Parent = this;
  // must be false to avoid reloading by vtkAppendCompositeDataLeaves::Update()
  this->Refresh = false;

  // INTIALIZE FILE NAME
  this->FileName = NULL;
  this->FileNameOld = new vtkStdString;

  // Case path
  this->CasePath = vtkCharArray::New();

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

  // for reading single precision binary format
  this->IsSinglePrecisionBinary = 0; // turned off by default
  this->IsSinglePrecisionBinaryOld = 0;

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
vtkNewOpenFOAMReader::~vtkNewOpenFOAMReader()
{
  this->LagrangianPaths->Delete();

  this->PatchDataArraySelection->Delete();
  this->CellDataArraySelection->Delete();
  this->PointDataArraySelection->Delete();
  this->LagrangianDataArraySelection->Delete();

  this->Readers->Delete();
  this->CasePath->Delete();

  this->SetFileName(0);
  delete this->FileNameOld;
}

//-----------------------------------------------------------------------------
// CanReadFile
int vtkNewOpenFOAMReader::CanReadFile(const char *vtkNotUsed(fileName))
{
  return 1; // so far CanReadFile does nothing.
}

//-----------------------------------------------------------------------------
// PrintSelf
void vtkNewOpenFOAMReader::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
  os << indent << "File Name: " << (this->FileName ? this->FileName : "(none)")
      << endl;
  os << indent << "Refresh: " << this->Refresh << endl;
  os << indent << "CreateCellToPoint: " << this->CreateCellToPoint << endl;
  os << indent << "CacheMesh: " << this->CacheMesh << endl;
  os << indent << "DecomposePolyhedra: " << this->DecomposePolyhedra << endl;
  os << indent << "PositionsIsIn13Format: " << this->PositionsIsIn13Format
      << endl;
  os << indent << "IsSinglePrecisionBinary: " << this->IsSinglePrecisionBinary
      << endl;
  os << indent << "ReadZones: " << this->ReadZones << endl;
  os << indent << "ListTimeStepsByControlDict: "
      << this->ListTimeStepsByControlDict << endl;
  os << indent << "AddDimensionsToArrayNames: "
      << this->AddDimensionsToArrayNames << endl;

  os << indent << "Case Path: \n";
  this->CasePath->PrintSelf(os, indent.GetNextIndent());

  this->Readers->InitTraversal();
  vtkObject *reader;
  while ((reader = this->Readers->GetNextItemAsObject()) != NULL)
    {
    os << indent << "Reader instance " << static_cast<void *>(reader) << ": \n";
    reader->PrintSelf(os, indent.GetNextIndent());
    }

  os << indent << "Patch Data Array Selection: \n";
  this->PatchDataArraySelection->PrintSelf(os, indent.GetNextIndent());
  os << indent << "Cell Data Array Selection: \n";
  this->CellDataArraySelection->PrintSelf(os, indent.GetNextIndent());
  os << indent << "Point Data Array Selection: \n";
  this->PointDataArraySelection->PrintSelf(os, indent.GetNextIndent());
  os << indent << "Lagrangian Data Array Selection: \n";
  this->LagrangianDataArraySelection->PrintSelf(os, indent.GetNextIndent());

  os << indent << "Patch Selection MTime Old: "
    << this->PatchSelectionMTimeOld << endl;
  os << indent << "Cell Selection MTime Old: "
    << this->CellSelectionMTimeOld << endl;
  os << indent << "Point Selection MTime Old: "
    << this->PointSelectionMTimeOld << endl;
  os << indent << "Lagrangian Selection MTime Old: "
    << this->LagrangianSelectionMTimeOld << endl;

  os << indent << "FileNameOld: " << (this->FileNameOld ? *this->FileNameOld :
      vtkStdString()) << endl;
  os << indent << "CreateCellToPointOld: " << this->CreateCellToPointOld
      << endl;
  os << indent << "DecomposePolyhedraOld: " << this->DecomposePolyhedraOld
      << endl;
  os << indent << "PositionsIsIn13FormatOld: " << this->PositionsIsIn13FormatOld
      << endl;
  os << indent << "IsSinglePrecisionBinaryOld: "
      << this->IsSinglePrecisionBinaryOld << endl;
  os << indent << "ReadZonesOld: " << this->ReadZonesOld << endl;
  os << indent << "ListTimeStepsByControlDictOld: "
      << this->ListTimeStepsByControlDictOld << endl;
  os << indent << "AddDimensionsToArrayNamesOld: "
      << this->AddDimensionsToArrayNamesOld << endl;

  os << indent << "Lagrangian Paths: \n";
  this->LagrangianPaths->PrintSelf(os, indent.GetNextIndent());

  os << indent << "Number Of Readers: " << this->NumberOfReaders << endl;
  os << indent << "Current Reader Index: " << this->CurrentReaderIndex << endl;

  os << indent << "Parent reader instance: "
      << static_cast<void *>(this->Parent) << endl;
  return;
}

//-----------------------------------------------------------------------------
// selection list handlers

int vtkNewOpenFOAMReader::GetNumberOfSelectionArrays(vtkDataArraySelection *s)
{
  return s->GetNumberOfArrays();
}

int vtkNewOpenFOAMReader::GetSelectionArrayStatus(vtkDataArraySelection *s,
    const char *name)
{
  return s->ArrayIsEnabled(name);
}

void vtkNewOpenFOAMReader::SetSelectionArrayStatus(vtkDataArraySelection *s,
    const char* name, int status)
{
  unsigned long int mTime = s->GetMTime();
  if (status)
    {
    s->EnableArray(name);
    }
  else
    {
    s->DisableArray(name);
    }
  if (mTime != s->GetMTime()) // indicate that the pipeline needs to be updated
    {
    this->Modified();
    }
}

const char *vtkNewOpenFOAMReader::GetSelectionArrayName(vtkDataArraySelection *s,
    int index)
{
  return s->GetArrayName(index);
}

void vtkNewOpenFOAMReader::DisableAllSelectionArrays(vtkDataArraySelection *s)
{
  unsigned long int mTime = s->GetMTime();
  s->DisableAllArrays();
  if (mTime != s->GetMTime())
    {
    this->Modified();
    }
}

void vtkNewOpenFOAMReader::EnableAllSelectionArrays(vtkDataArraySelection *s)
{
  unsigned long int mTime = s->GetMTime();
  s->EnableAllArrays();
  if (mTime != s->GetMTime())
    {
    this->Modified();
    }
}

//-----------------------------------------------------------------------------
// RequestInformation
int vtkNewOpenFOAMReader::RequestInformation(vtkInformation *vtkNotUsed(request), vtkInformationVector **vtkNotUsed(inputVector), vtkInformationVector *outputVector)
{
  if (!this->FileName || strlen(this->FileName) == 0)
    {
    vtkErrorMacro("FileName has to be specified!");
    return 0;
    }

  if (this->Parent == this && (*this->FileNameOld != this->FileName
      || this->ListTimeStepsByControlDict
          != this->ListTimeStepsByControlDictOld || this->Refresh))
    {
    // retain selection status when just refreshing a case
    if (*this->FileNameOld != "" && *this->FileNameOld != this->FileName)
      {
      // clear selections
      this->CellDataArraySelection->RemoveAllArrays();
      this->PointDataArraySelection->RemoveAllArrays();
      this->LagrangianDataArraySelection->RemoveAllArrays();
      this->PatchDataArraySelection->RemoveAllArrays();
      }

    // Reset NumberOfReaders here so that the variable will not be
    // reset unwantedly when MakeInformationVector() is called from
    // vtkNewPOpenFOAMReader
    this->NumberOfReaders = 0;

    if (!this->MakeInformationVector(outputVector, vtkStdString(""))
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
int vtkNewOpenFOAMReader::RequestData(vtkInformation *vtkNotUsed(request), vtkInformationVector **vtkNotUsed(inputVector), vtkInformationVector *outputVector)
{
  vtkInformation* outInfo = outputVector->GetInformationObject(0);
  vtkMultiBlockDataSet
      *output =
          vtkMultiBlockDataSet::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

  int nSteps = 0;
  double *requestedTimeValues = NULL;
  if (outInfo->Has(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEPS()))
    {
    requestedTimeValues
        = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEPS());
    nSteps = outInfo->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
    }

  if (nSteps > 0)
    {
    outInfo->Set(vtkDataObject::DATA_TIME_STEPS(), requestedTimeValues, 1);
    this->SetTimeValue(requestedTimeValues[0]);
    }

  if (this->Parent == this)
    {
    output->GetFieldData()->AddArray(this->CasePath);
    if (!this->MakeMetaDataAtTimeStep(false))
      {
      return 0;
      }
    this->CurrentReaderIndex = 0;
    }

  // compute flags
  // internal mesh selection change is detected within each reader
  const bool recreateInternalMesh = (!this->Parent->CacheMesh)
      || this->Parent->DecomposePolyhedra
          != this->Parent->DecomposePolyhedraOld || this->Parent->ReadZones
      != this->Parent->ReadZonesOld || this->Parent->ListTimeStepsByControlDict
      != this->Parent->ListTimeStepsByControlDictOld
      || this->Parent->IsSinglePrecisionBinary
          != this->Parent->IsSinglePrecisionBinaryOld;
  const bool recreateBoundaryMesh =
      this->Parent->PatchDataArraySelection->GetMTime()
          != this->Parent->PatchSelectionMTimeOld
          || this->Parent->CreateCellToPoint
              != this->Parent->CreateCellToPointOld;
  const bool updateVariables = this->Parent->CellDataArraySelection->GetMTime()
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
  vtkNewOpenFOAMReaderPrivate *reader;
  // if the only region is not a subregion, omit being wrapped by a
  // multiblock dataset
  if (this->Readers->GetNumberOfItems() == 1 && (reader = vtkNewOpenFOAMReaderPrivate::SafeDownCast(
          this->Readers->GetItemAsObject(0)))->GetRegionName() == "")
    {
    ret = reader->RequestData(output, recreateInternalMesh,
        recreateBoundaryMesh, updateVariables);
    this->Parent->CurrentReaderIndex++;
    }
  else
    {
    this->Readers->InitTraversal();
    while ((reader
        = vtkNewOpenFOAMReaderPrivate::SafeDownCast(this->Readers->GetNextItemAsObject()))
        != NULL)
      {
      vtkMultiBlockDataSet *subOutput = vtkMultiBlockDataSet::New();
      if (reader->RequestData(subOutput, recreateInternalMesh,
          recreateBoundaryMesh, updateVariables))
        {
        vtkStdString regionName(reader->GetRegionName());
        if (regionName == "")
          {
          regionName = "defaultRegion";
          }
        const int blockI = output->GetNumberOfBlocks();
        output->SetBlock(blockI, subOutput);
        output->GetMetaData(blockI)->Set(vtkCompositeDataSet::NAME(), regionName.c_str());
        }
      else
        {
        ret = 0;
        }
      subOutput->Delete();
      this->Parent->CurrentReaderIndex++;
      }
    }

  if (this->Parent == this) // update only if this is the top-level reader
    {
    this->UpdateStatus();
    }

  return ret;
}

//-----------------------------------------------------------------------------
void vtkNewOpenFOAMReader::SetTimeInformation(vtkInformationVector *outputVector,
    vtkDoubleArray *timeValues)
{
  if (timeValues->GetNumberOfTuples() > 0)
    {
    outputVector->GetInformationObject(0)->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(),
        timeValues->GetPointer(0), timeValues->GetNumberOfTuples());

    double timeRange[2];
    timeRange[0] = timeValues->GetValue(0);
    timeRange[1] = timeValues->GetValue(timeValues->GetNumberOfTuples() - 1);
    outputVector->GetInformationObject(0)->Set(vtkStreamingDemandDrivenPipeline::TIME_RANGE(), timeRange, 2);
    }
  else
    {
    double timeRange[2];
    timeRange[0] = timeRange[1] = 0.0;
    outputVector->GetInformationObject(0)->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(), timeRange, 0);
    outputVector->GetInformationObject(0)->Set(vtkStreamingDemandDrivenPipeline::TIME_RANGE(), timeRange, 2);
    }
}

//-----------------------------------------------------------------------------
void vtkNewOpenFOAMReader::GetRegions(vtkStringArray *regionNames,
    const vtkStdString &timeDir)
{
  vtkDirectory *dir = vtkDirectory::New();
  if (!dir->Open(timeDir.c_str()))
    {
    dir->Delete();
    return;
    }
  for (int fileI = 0; fileI < dir->GetNumberOfFiles(); fileI++)
    {
    vtkStdString subDir(dir->GetFile(fileI));
    // "uniform" subdir is used for recording time by OF <= 1.4.1
    if (dir->FileIsDirectory(subDir.c_str()) && subDir != "." && subDir != ".."
        && subDir != "polyMesh" && subDir != "lagrangian" && subDir != "uniform")
      {
      if(regionNames->LookupValue(subDir) >= 0)
        {
        continue;
        }
      vtkStdString boundaryPath(timeDir + "/" + subDir + "/polyMesh/boundary");
      if (vtksys::SystemTools::FileExists(boundaryPath.c_str(), true)
          || vtksys::SystemTools::FileExists((boundaryPath + ".gz").c_str(), true))
        {
        regionNames->InsertNextValue(subDir);
        }
      }
    }
  dir->Delete();
}

//-----------------------------------------------------------------------------
int vtkNewOpenFOAMReader::MakeInformationVector(
    vtkInformationVector *outputVector, const vtkStdString& procName)
{
  *this->FileNameOld = vtkStdString(this->FileName);

  // clear prior case information
  this->Readers->RemoveAllItems();

  // recreate case information
  vtkStdString casePath, controlDictPath;
  this->CreateCasePath(casePath, controlDictPath);
  casePath += procName + (procName == "" ? "" : "/");
  vtkNewOpenFOAMReaderPrivate *masterReader = vtkNewOpenFOAMReaderPrivate::New();
  if (!masterReader->MakeInformationVector(casePath, controlDictPath, procName,
      this->Parent))
    {
    masterReader->Delete();
    return 0;
    }

  if (masterReader->GetTimeValues()->GetNumberOfTuples() == 0)
    {
    vtkErrorMacro(<< "Case " << casePath.c_str()
        << " contains no timestep data.");
    masterReader->Delete();
    return 0;
    }

  this->Readers->AddItem(masterReader);

  if (outputVector != NULL)
    {
    this->SetTimeInformation(outputVector, masterReader->GetTimeValues());
    }

  // search subregions under constant subdirectory
  vtkStringArray *regionNames = vtkStringArray::New();
  this->GetRegions(regionNames, casePath + "constant");
  vtkStringArray *timeNames = masterReader->GetTimeNames();
  int nTimes = timeNames->GetNumberOfValues() > 2
      ? 2 : timeNames->GetNumberOfValues();
  for(int timeI = 0; timeI < nTimes; timeI++)
    {
    this->GetRegions(regionNames, casePath + timeNames->GetValue(timeI));
    }
  regionNames->Squeeze();
  for(int regionI = 0; regionI < regionNames->GetNumberOfValues(); regionI++)
    {
    vtkNewOpenFOAMReaderPrivate *subReader = vtkNewOpenFOAMReaderPrivate::New();
    subReader->SetupInformation(casePath, regionNames->GetValue(regionI),
        procName, masterReader);
    this->Readers->AddItem(subReader);
    subReader->Delete();
    }
  regionNames->Delete();
  masterReader->Delete();
  this->Parent->NumberOfReaders += this->Readers->GetNumberOfItems();

  if (this->Parent == this)
    {
    this->CreateCharArrayFromString(this->CasePath, "CasePath", casePath);
    }

  return 1;
}

//-----------------------------------------------------------------------------
void vtkNewOpenFOAMReader::CreateCasePath(vtkStdString &casePath,
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
  if (pos == vtkStdString::npos)
    {
    // if there's no prepending path, prefix with the current directory
    controlDictPath = "." + pathSeparator + controlDictPath;
    pos = 1;
    }
  if (controlDictPath.substr(pos + 1, 11) == "controlDict")
    {
    // remove trailing "/controlDict*"
    casePath = controlDictPath.substr(0, pos - 1);
    if (casePath == ".")
      {
      casePath = ".." + pathSeparator;
      }
    else
      {
      pos = casePath.find_last_of(pathFindSeparator);
      if (pos == vtkStdString::npos)
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
void vtkNewOpenFOAMReader::AddSelectionNames(vtkDataArraySelection *selections,
    vtkStringArray *objects)
{
  objects->Squeeze();
  vtkSortDataArray::Sort(objects);
  for (int nameI = 0; nameI < objects->GetNumberOfValues(); nameI++)
    {
    selections->AddArray(objects->GetValue(nameI).c_str());
    }
  objects->Delete();
}

//-----------------------------------------------------------------------------
bool vtkNewOpenFOAMReader::SetTimeValue(const double timeValue)
{
  bool modified = false;
  vtkNewOpenFOAMReaderPrivate *reader;
  this->Readers->InitTraversal();
  while ((reader
      = vtkNewOpenFOAMReaderPrivate::SafeDownCast(this->Readers->GetNextItemAsObject()))
      != NULL)
    {
    const unsigned long mTime = reader->GetMTime();
    reader->SetTimeValue(timeValue);
    if (reader->GetMTime() != mTime)
      {
      modified = true;
      }
    }
  return modified;
}

//-----------------------------------------------------------------------------
vtkDoubleArray *vtkNewOpenFOAMReader::GetTimeValues()
{
  if (this->Readers->GetNumberOfItems() <= 0)
    {
    return NULL;
    }
  vtkNewOpenFOAMReaderPrivate *reader =
      vtkNewOpenFOAMReaderPrivate::SafeDownCast(this->Readers->GetItemAsObject(0));
  return reader != NULL ? reader->GetTimeValues() : NULL;
}

//-----------------------------------------------------------------------------
int vtkNewOpenFOAMReader::MakeMetaDataAtTimeStep(const bool listNextTimeStep)
{
  vtkStringArray *cellSelectionNames = vtkStringArray::New();
  vtkStringArray *pointSelectionNames = vtkStringArray::New();
  vtkStringArray *lagrangianSelectionNames = vtkStringArray::New();
  int ret = 1;
  vtkNewOpenFOAMReaderPrivate *reader;
  this->Readers->InitTraversal();
  while ((reader
      = vtkNewOpenFOAMReaderPrivate::SafeDownCast(this->Readers->GetNextItemAsObject()))
      != NULL)
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
void vtkNewOpenFOAMReader::CreateCharArrayFromString(vtkCharArray *array,
    const char *name, vtkStdString &string)
{
  array->Initialize();
  array->SetName(name);
  const size_t len = string.length();
  char *ptr = array->WritePointer(0, static_cast<vtkIdType>(len + 1));
  memcpy(ptr, string.c_str(), len);
  ptr[len] = '\0';
}

//-----------------------------------------------------------------------------
void vtkNewOpenFOAMReader::UpdateStatus()
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
  this->IsSinglePrecisionBinaryOld = this->IsSinglePrecisionBinary;
  this->ReadZonesOld = this->ReadZones;
  this->ListTimeStepsByControlDictOld = this->ListTimeStepsByControlDict;
  this->AddDimensionsToArrayNamesOld = this->AddDimensionsToArrayNames;
}

//-----------------------------------------------------------------------------
void vtkNewOpenFOAMReader::UpdateProgress(double amount)
{
  this->vtkAlgorithm::UpdateProgress((static_cast<double>(this->Parent->CurrentReaderIndex)
      + amount) / static_cast<double>(this->Parent->NumberOfReaders));
}

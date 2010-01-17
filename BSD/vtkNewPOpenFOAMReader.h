/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkNewPOpenFOAMReader.h,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See License_v1.2.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkNewPOpenFOAMReader - reads a decomposed dataset in OpenFOAM(R) format
// .SECTION Description
// vtkNewPOpenFOAMReader creates a multiblock dataset. It reads
// parallel-decomposed mesh information and time dependent data.  The
// polyMesh folders contain mesh information. The time folders contain
// transient data for the cells. Each folder can contain any number of
// data files.

// .SECTION Thanks
// This class was developed by Takuya Oshima at Niigata University,
// Japan (oshima@eng.niigata-u.ac.jp).
// OPENFOAM(R) is a registered trade mark of OpenCFD Limited, the
// producer of the OpenFOAM software and owner of the OPENFOAM(R) and
// OpenCFD(R) trade marks. This code is not approved or endorsed by
// OpenCFD Limited.



#ifndef __vtkNewPOpenFOAMReader_h
#define __vtkNewPOpenFOAMReader_h

#include "vtkNewOpenFOAMReader.h"

class vtkDataArraySelection;
class vtkMultiProcessController;

class
#if !defined(POpenFOAMReaderPlugin_EXPORTS)
VTK_PARALLEL_EXPORT
#endif
vtkNewPOpenFOAMReader : public vtkNewOpenFOAMReader
{
public:
  //BTX
  enum caseType { DECOMPOSED_CASE = 0, RECONSTRUCTED_CASE = 1 };
  //ETX
  static vtkNewPOpenFOAMReader *New();

  vtkTypeRevisionMacro(vtkNewPOpenFOAMReader, vtkNewOpenFOAMReader);

  void PrintSelf(ostream &os, vtkIndent indent);

  // Description:
  // Set and get case type. 0 = decomposed case, 1 = reconstructed case.
  void SetCaseType(const int t);
  vtkGetMacro(CaseType, caseType);
  // Description:
  // Set and get the controller.
  virtual void SetController(vtkMultiProcessController *);
  vtkGetObjectMacro(Controller, vtkMultiProcessController);

  // Description:
  // Set and get the polling interval for watching the case directory.
  // Not using vtkSetMacro so as not to modify MTime.
  virtual void SetUiInterval(const int interval)
    { this->UiInterval = interval; }
  vtkGetMacro(UiInterval, int);

  // Description:
  // Set and get if rescale the range of color lookup table
  // automatically on refresh.
  // Not using vtkSetMacro so as not to modify MTime.
  virtual void SetUiRescale(const int rescale) { this->UiRescale = rescale; }
  vtkGetMacro(UiRescale, int);
  vtkBooleanMacro(UiRescale, int);

  // Description:
  // Set and get if rescale the range of color lookup table
  // automatically on refresh.
  // Not using vtkSetMacro so as not to modify MTime.
  virtual void SetUiWatch(const int watch) { this->UiWatch = watch; }
  vtkGetMacro(UiWatch, int);
  vtkBooleanMacro(UiWatch, int);

protected:
  vtkNewPOpenFOAMReader();
  ~vtkNewPOpenFOAMReader();

  int RequestInformation(vtkInformation *, vtkInformationVector **,
    vtkInformationVector *);
  int RequestData(vtkInformation *, vtkInformationVector **,
    vtkInformationVector *);

private:
  vtkMultiProcessController *Controller;
  caseType CaseType;
  unsigned long MTimeOld;
  int MaximumNumberOfPieces;
  int NumProcesses;
  int ProcessId;

  // the followings are for additional UI elements
  int UiInterval;
  int UiRescale;
  int UiWatch;

  vtkNewPOpenFOAMReader(const vtkNewPOpenFOAMReader &); // Not implemented.
  void operator=(const vtkNewPOpenFOAMReader &); // Not implemented.

  void GatherMetaData();
  void BroadcastStatus(int &);
  void Broadcast(vtkStringArray *);
  int ConstructBlocks(vtkMultiBlockDataSet *, int *, int);
  void BroadcastStructure(vtkMultiBlockDataSet *, const int);
  void AllGather(vtkStringArray *);
  void AllGather(vtkDataArraySelection *);
};

#endif

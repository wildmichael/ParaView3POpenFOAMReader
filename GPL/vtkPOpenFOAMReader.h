/*=========================================================================

    This file is part of vtkPOpenFOAMReader.

    vtkPOpenFOAMReader is free software: you can redistribute it
    and/or modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation, either version 3 of
    the License, or (at your option) any later version.

    vtkPOpenFOAMReader is distributed in the hope that it will be
    useful, but WITHOUT ANY WARRANTY; without even the implied
    warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
    See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with vtkPOpenFOAMReader.  If not, see
    <http://www.gnu.org/licenses/>.

=========================================================================*/
// Copyright (c) 2008 Takuya OSHIMA <oshima@eng.niigata-u.ac.jp>.
// All rights reserved.
// Date: 2008-09-27

// .NAME vtkPOpenFOAMReader - reads a decomposed dataset in OpenFOAM format
// .SECTION Description
// vtkPOpenFOAMReader creates a multiblock dataset. It reads
// parallel-decomposed mesh information and time dependent data.  The
// polyMesh folders contain mesh information. The time folders contain
// transient data for the cells. Each folder can contain any number of
// data files.

// .SECTION Thanks
// This class was developed by Takuya Oshima at Niigata University,
// Japan (oshima@eng.niigata-u.ac.jp).

#ifndef __vtkPOpenFOAMReader_h
#define __vtkPOpenFOAMReader_h

#include "vtkOpenFOAMReader.h"

class vtkDataArraySelection;
class vtkMultiProcessController;

class VTK_EXPORT vtkPOpenFOAMReader : public vtkOpenFOAMReader
{
public:
  //BTX
  enum caseType { DECOMPOSED_CASE = 0, RECONSTRUCTED_CASE = 1 };
  //ETX
  static vtkPOpenFOAMReader *New();
  vtkTypeRevisionMacro(vtkPOpenFOAMReader, vtkOpenFOAMReader);

  void PrintSelf(ostream &os, vtkIndent indent);

  // Description:
  // Set and get case type. 0 = decomposed case, 1 = reconstructed case.
  void SetCaseType(const int t);
  vtkGetMacro(CaseType, caseType);
  // Description:
  // Set and get the controller.
  virtual void SetController(vtkMultiProcessController *);
  vtkGetObjectMacro(Controller, vtkMultiProcessController);

protected:
  vtkPOpenFOAMReader();
  ~vtkPOpenFOAMReader();

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

  vtkPOpenFOAMReader(const vtkPOpenFOAMReader &); // Not implemented.
  void operator=(const vtkPOpenFOAMReader &); // Not implemented.

  void GatherMetaData();
  void BroadcastStatus(int &);
  void Broadcast(vtkStringArray *);
  //void Gather(vtkDataArraySelection *);
  void AllGather(vtkStringArray *);
  void AllGather(vtkDataArraySelection *);
};

#endif

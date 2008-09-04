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
// Date: 2008-08-31

#include "vtkPOpenFOAMReader.h"

#include "vtkAppendCompositeDataLeaves.h"
#include "vtkCollection.h"
#include "vtkDataArraySelection.h"
#include "vtkDirectory.h"
#include "vtkDoubleArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkIntArray.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkMultiProcessController.h"
#include "vtkObjectFactory.h"
#include "vtkSortDataArray.h"
#include "vtkStdString.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkStringArray.h"

vtkCxxRevisionMacro(vtkPOpenFOAMReader, "$Revision: 1.00 $");
vtkStandardNewMacro(vtkPOpenFOAMReader);

//-----------------------------------------------------------------------------
vtkPOpenFOAMReader::vtkPOpenFOAMReader()
{
  this->CaseType = DECOMPOSED_CASE;
  this->MTimeOld = 0;
  this->MaximumNumberOfPieces = 1;
}

//-----------------------------------------------------------------------------
vtkPOpenFOAMReader::~vtkPOpenFOAMReader()
{
}

//-----------------------------------------------------------------------------
void vtkPOpenFOAMReader::PrintSelf(ostream &os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
  os << indent << "MTimeOld: " << this->MTimeOld;
  os << indent << "MaximumNumberOfPieces: " << this->MaximumNumberOfPieces;
}

//-----------------------------------------------------------------------------
void vtkPOpenFOAMReader::SetCaseType(const int t)
{
  if(this->CaseType != t)
    {
    this->CaseType = static_cast<caseType>(t);
    this->Refresh = true;
    this->Modified();
    }
}

//-----------------------------------------------------------------------------
int vtkPOpenFOAMReader::RequestInformation(vtkInformation *request,
  vtkInformationVector **inputVector, vtkInformationVector *outputVector)
{
  if(this->CaseType == RECONSTRUCTED_CASE)
    {
    return this->Superclass::RequestInformation(request, inputVector,
      outputVector);
    }

  if(!this->Superclass::FileName || strlen(this->Superclass::FileName) == 0)
    {
    vtkErrorMacro("FileName has to be specified!");
    return 0;
    }

  if(*this->Superclass::FileNameOld != vtkStdString(this->Superclass::FileName)
    || this->Superclass::ListTimeStepsByControlDict
    != this->Superclass::ListTimeStepsByControlDictOld
    || this->Superclass::Refresh)
    {
    *this->Superclass::FileNameOld = vtkStdString(this->FileName);

    // clear prior case information
    this->Superclass::Readers->RemoveAllItems();
    this->Superclass::NumberOfReaders = 0;

    this->Superclass::CellDataArraySelection->RemoveAllArrays();
    this->Superclass::PointDataArraySelection->RemoveAllArrays();
    this->Superclass::LagrangianDataArraySelection->RemoveAllArrays();
    this->Superclass::PatchDataArraySelection->RemoveAllArrays();

    vtkMultiProcessController *ctrl
      = vtkMultiProcessController::GetGlobalController();
    const int processId = (ctrl == NULL ? 0 : ctrl->GetLocalProcessId());

    vtkStringArray *procNames = vtkStringArray::New();
    vtkDoubleArray *timeValues;
    int ret = 1;
    if(processId == 0)
      {
      // recreate case information
      vtkStdString masterCasePath, controlDictPath;
      this->Superclass::CreateCasePath(masterCasePath, controlDictPath);

      // search and list processor subdirectories
      vtkDirectory *dir = vtkDirectory::New();
      if(!dir->Open(masterCasePath.c_str()))
        {
        vtkErrorMacro(<< "Can't open " << masterCasePath.c_str());
        dir->Delete();
        ctrl->Broadcast(&(ret = 0), 1, 0);
        return 0;
        }
      vtkIntArray *procNos = vtkIntArray::New();
      for(int fileI = 0; fileI < dir->GetNumberOfFiles(); fileI++)
        {
        const vtkStdString subDir(dir->GetFile(fileI));
        if(subDir.substr(0, 9) == "processor")
          {
          const vtkStdString procNoStr(subDir.substr(9, vtkStdString::npos));
          char *conversionEnd;
          const int procNo = strtol(procNoStr.c_str(), &conversionEnd, 10);
          if(procNoStr.c_str() + procNoStr.length() == conversionEnd
            && procNo >= 0)
            {
            procNos->InsertNextValue(procNo);
            procNames->InsertNextValue(subDir);
            }
          }
        }
      procNos->Squeeze();
      procNames->Squeeze();
      dir->Delete();

#if 0
      if(procNos->GetNumberOfTuples() == 0)
        {
        procNos->Delete();
        procNames->Delete();
        ctrl->Broadcast(&(ret = 0), 1, 0);
        vtkErrorMacro(<< "No processor subdirectory found.");
        return 0;
        }
#endif
      // sort processor subdirectories by processor numbers
      vtkSortDataArray::Sort(procNos, procNames);
      procNos->Delete();

      // get time directories from the first processor subdirectory
      if(procNames->GetNumberOfTuples() > 0)
	{
	vtkOpenFOAMReader *masterReader = vtkOpenFOAMReader::New();
        masterReader->SetFileName(this->FileName);
        masterReader->SetParent(this);
        if(!masterReader->MakeInformationVector(outputVector, procNames
          ->GetValue(0)) || !masterReader->MakeMetaDataAtTimeStep(true))
          {
          procNames->Delete();
          masterReader->Delete();
          ctrl->Broadcast(&(ret = 0), 1, 0);
          return 0;
          }
        this->Superclass::Readers->AddItem(masterReader);
        timeValues = masterReader->GetTimeValues();
        masterReader->Delete();
	}
      else
	{
        timeValues = vtkDoubleArray::New();
	this->Superclass::SetTimeInformation(outputVector, timeValues);
	}
      }
    else
      {
      timeValues = vtkDoubleArray::New();
      }

    const int nProcesses = (ctrl == NULL ? 1 : ctrl->GetNumberOfProcesses());
    if(nProcesses > 1)
      {
      // if there was an error in process 0 abort all processes
      ctrl->Broadcast(&ret, 1, 0);
      if(ret == 0)
        {
        vtkErrorMacro(<< "The master process returned an error.");
        timeValues->Delete(); // don't have to care about process 0
        return 0;
        }

      this->Broadcast(procNames, ctrl);
      ctrl->Broadcast(timeValues, 0);
      if(processId != 0)
        {
        this->Superclass::SetTimeInformation(outputVector, timeValues);
        timeValues->Delete();
        }
      }

    if(processId == 0 && procNames->GetNumberOfTuples() == 0)
      {
      timeValues->Delete();
      }

    this->MaximumNumberOfPieces = procNames->GetNumberOfTuples();

    // create reader instances for other processor subdirectories
    for(int procI = (processId == 0 ? nProcesses : processId);
      procI < procNames->GetNumberOfTuples(); procI += nProcesses)
      {
      vtkOpenFOAMReader *subReader = vtkOpenFOAMReader::New();
      subReader->SetFileName(this->FileName);
      subReader->SetParent(this);
      // if getting metadata failed simply delete the reader instance
      if(subReader->MakeInformationVector(NULL, procNames->GetValue(procI))
        && subReader->MakeMetaDataAtTimeStep(true))
        {
        this->Superclass::Readers->AddItem(subReader);
        }
      else
        {
        vtkWarningMacro(<<"Removing reader for processor subdirectory "
          << procNames->GetValue(procI).c_str());
        }
      subReader->Delete();
      }

    procNames->Delete();

    this->GatherMetaData();
    this->Superclass::Refresh = false;
    }

  // it seems MAXIMUM_NUMBER_OF_PIECES must be returned every time
  // RequestInformation() is called
  outputVector->GetInformationObject(0)->Set(
    vtkStreamingDemandDrivenPipeline::MAXIMUM_NUMBER_OF_PIECES(),
    this->MaximumNumberOfPieces);

  return 1;
}

//-----------------------------------------------------------------------------
int vtkPOpenFOAMReader::RequestData(vtkInformation *request,
  vtkInformationVector **inputVector, vtkInformationVector *outputVector)
{
  if(this->CaseType == RECONSTRUCTED_CASE)
    {
    return this->Superclass::RequestData(request, inputVector, outputVector);
    }

  vtkInformation* outInfo = outputVector->GetInformationObject(0);
  vtkMultiBlockDataSet *output = vtkMultiBlockDataSet::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  int ret = 1;
  if(this->Superclass::Readers->GetNumberOfItems() > 0)
    {
    int nSteps = 0;
    double *requestedTimeValues = NULL;
    if(outInfo->Has(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEPS()))
      {
      requestedTimeValues
        = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEPS());
      nSteps = outInfo->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
      if(nSteps > 0)
        {
        outInfo->Set(vtkDataObject::DATA_TIME_STEPS(), requestedTimeValues, 1);
        }
      }
#if 0
    const int nPieces = outInfo->Get(
      vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES());
    // returns -1 for processes whose processId exceeds MAXIMUM_NUMBER_OF_PIECES
    // which had been set in RequestInformation()
    const int nMaxPieces = outInfo->Get(
      vtkStreamingDemandDrivenPipeline::MAXIMUM_NUMBER_OF_PIECES());
    const int piece
      = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
#endif
    vtkAppendCompositeDataLeaves *append = vtkAppendCompositeDataLeaves::New();
    // append->AppendFieldDataOn();

    vtkOpenFOAMReader *reader;
    this->Superclass::CurrentReaderIndex = 0;
    this->Superclass::Readers->InitTraversal();
    while((reader = vtkOpenFOAMReader::SafeDownCast(
      this->Superclass::Readers->GetNextItemAsObject())) != NULL)
      {
      // even if the child readers themselves are not modified, mark
      // them as modified if "this" has been modified, since they
      // refer to the property of "this"
      if((nSteps > 0 && reader->SetTimeValue(requestedTimeValues[0]))
        || this->MTimeOld != this->GetMTime())
        {
        reader->Modified();
        }
      if(reader->MakeMetaDataAtTimeStep(false))
        {
        append->AddInputConnection(reader->GetOutputPort());
        }
      }

    this->GatherMetaData();

    if(append->GetInput() == NULL)
      {
      output->Initialize();
      ret = 0;
      }
    else
      {
      // reader->RequestInformation() and RequestData() are called
      // for all reader instances without setting UPDATE_TIME_STEPS
      append->Update();
      output->ShallowCopy(append->GetOutput());
      }
    append->Delete();
    }
  else
    {
    this->GatherMetaData();
    // page 322 of The ParaView Guide says the output must be initialized
    output->Initialize();
    }

  this->Superclass::UpdateStatus();
  this->MTimeOld = this->GetMTime();

  return ret;
}

//-----------------------------------------------------------------------------
void vtkPOpenFOAMReader::GatherMetaData()
{
  vtkMultiProcessController *ctrl
    = vtkMultiProcessController::GetGlobalController();
  if(ctrl->GetNumberOfProcesses() > 1)
    {
    this->Gather(this->Superclass::PatchDataArraySelection, ctrl);
    this->Gather(this->Superclass::CellDataArraySelection, ctrl);
    this->Gather(this->Superclass::PointDataArraySelection, ctrl);
    this->Gather(this->Superclass::LagrangianDataArraySelection, ctrl);
    // omit removing duplicated entries of LagrangianPaths as well
    // when the number of processes is 1 assuming there's no duplicate
    // entry within a process
    this->AllGather(this->Superclass::LagrangianPaths, ctrl);
    }
}

//-----------------------------------------------------------------------------
// Broadcast a vtkStringArray on process 0 to all processes
void vtkPOpenFOAMReader::Broadcast(vtkStringArray *sa,
  vtkMultiProcessController *ctrl)
{
  vtkIdType lengths[2];
  if(ctrl->GetLocalProcessId() == 0)
    {
    lengths[0] = sa->GetNumberOfTuples();
    lengths[1] = 0;
    for(int strI = 0; strI < sa->GetNumberOfTuples(); strI++)
      {
      lengths[1] += sa->GetValue(strI).length() + 1;
      }
    }
  ctrl->Broadcast(lengths, 2, 0);
  char *contents = new char [lengths[1]];
  if(ctrl->GetLocalProcessId() == 0)
    {
    for(int strI = 0, idx = 0; strI < sa->GetNumberOfTuples(); strI++)
      {
      const int len = sa->GetValue(strI).length() + 1;
      memmove(contents + idx, sa->GetValue(strI).c_str(), len);
      idx += len;
      }
    }
  ctrl->Broadcast(contents, lengths[1], 0);
  if(ctrl->GetLocalProcessId() != 0)
    {
    sa->Initialize();
    sa->SetNumberOfTuples(lengths[0]);
    for(int strI = 0, idx = 0; strI < lengths[0]; strI++)
      {
      sa->SetValue(strI, contents + idx);
      idx += sa->GetValue(strI).length() + 1;
      }
    }
  delete [] contents;
}

//-----------------------------------------------------------------------------
// Gather vtkDataArraySelections on all processes to process 0
void vtkPOpenFOAMReader::Gather(vtkDataArraySelection *s,
  vtkMultiProcessController *ctrl)
{
  vtkIdType length = 0;
  for(int strI = 0; strI < s->GetNumberOfArrays(); strI++)
    {
    length += strlen(s->GetArrayName(strI)) + 2;
    }
  vtkIdType *lengths = NULL;
  if(ctrl->GetLocalProcessId() == 0)
    {
    lengths = new vtkIdType [ctrl->GetNumberOfProcesses()];
    }
  ctrl->Gather(&length, lengths, 1, 0);
  char *allContents = NULL;
  vtkIdType totalLength = 0, *offsets = NULL;
  if(ctrl->GetLocalProcessId() == 0)
    {
    offsets = new vtkIdType [ctrl->GetNumberOfProcesses()];
    for(int procI = 0; procI < ctrl->GetNumberOfProcesses(); procI++)
      {
      offsets[procI] = totalLength;
      totalLength += lengths[procI];
      }
    allContents = new char [totalLength];
    }
  char *contents = new char [length];
  for(int strI = 0, idx = 0; strI < s->GetNumberOfArrays(); strI++)
    {
    const char *arrayName = s->GetArrayName(strI);
    contents[idx] = s->ArrayIsEnabled(arrayName);
    const int len = strlen(arrayName) + 1;
    memmove(contents + idx + 1, arrayName, len);
    idx += len + 1;
    }
  ctrl->GatherV(contents, allContents, length, lengths, offsets, 0);
  delete [] contents;
  if(ctrl->GetLocalProcessId() == 0)
    {
    delete [] lengths;
    delete [] offsets;
    // s->RemoveAllArrays();
    for(int idx = 0; idx < totalLength;
      idx += strlen(allContents + idx + 1) + 2)
      {
      const char *arrayName = allContents + idx + 1;
      // set array status only when the array is not existent
      if(!s->ArrayExists(arrayName))
        {
        if(allContents[idx] == 0)
          {
          s->DisableArray(arrayName);
          }
        else
          {
          s->EnableArray(arrayName);
          }
        }
      }
    delete [] allContents;
    }
}

//-----------------------------------------------------------------------------
// AllGather vtkStringArray from and to all processes
void vtkPOpenFOAMReader::AllGather(vtkStringArray *s,
  vtkMultiProcessController *ctrl)
{
  vtkIdType length = 0;
  for(int strI = 0; strI < s->GetNumberOfTuples(); strI++)
    {
    length += s->GetValue(strI).length() + 1;
    }
  vtkIdType *lengths = new vtkIdType [ctrl->GetNumberOfProcesses()];
  ctrl->AllGather(&length, lengths, 1);
  vtkIdType totalLength = 0;
  vtkIdType *offsets = new vtkIdType [ctrl->GetNumberOfProcesses()];
  for(int procI = 0; procI < ctrl->GetNumberOfProcesses(); procI++)
    {
    offsets[procI] = totalLength;
    totalLength += lengths[procI];
    }
  char *allContents = new char [totalLength], *contents = new char [length];
  for(int strI = 0, idx = 0; strI < s->GetNumberOfTuples(); strI++)
    {
    const int len = s->GetValue(strI).length() + 1;
    memmove(contents + idx, s->GetValue(strI).c_str(), len);
    idx += len;
    }
  ctrl->AllGatherV(contents, allContents, length, lengths, offsets);
  delete [] contents;
  delete [] lengths;
  delete [] offsets;
  s->Initialize();
  for(int idx = 0; idx < totalLength; idx += strlen(allContents + idx) + 1)
    {
    const char *str = allContents + idx;
    // insert only when the same string is not found
    if(s->LookupValue(str) == -1)
      {
      s->InsertNextValue(str);
      }
    }
  s->Squeeze();
  delete [] allContents;
}

//-----------------------------------------------------------------------------
// AllGather vtkDataArraySelections from and to all processes
void vtkPOpenFOAMReader::AllGather(vtkDataArraySelection *s,
  vtkMultiProcessController *ctrl)
{
  vtkIdType length = 0;
  for(int strI = 0; strI < s->GetNumberOfArrays(); strI++)
    {
    length += strlen(s->GetArrayName(strI)) + 2;
    }
  vtkIdType *lengths = new vtkIdType [ctrl->GetNumberOfProcesses()];
  ctrl->AllGather(&length, lengths, 1);
  vtkIdType totalLength = 0;
  vtkIdType *offsets = new vtkIdType [ctrl->GetNumberOfProcesses()];
  for(int procI = 0; procI < ctrl->GetNumberOfProcesses(); procI++)
    {
    offsets[procI] = totalLength;
    totalLength += lengths[procI];
    }
  char *allContents = new char [totalLength], *contents = new char [length];
  for(int strI = 0, idx = 0; strI < s->GetNumberOfArrays(); strI++)
    {
    const char *arrayName = s->GetArrayName(strI);
    contents[idx] = s->ArrayIsEnabled(arrayName);
    const int len = strlen(arrayName) + 1;
    memmove(contents + idx + 1, arrayName, len);
    idx += len + 1;
    }
  ctrl->AllGatherV(contents, allContents, length, lengths, offsets);
  delete [] contents;
  delete [] lengths;
  delete [] offsets;
  // s->RemoveAllArrays();
  for(int idx = 0; idx < totalLength; idx += strlen(allContents + idx + 1) + 2)
    {
    const char *arrayName = allContents + idx + 1;
    s->AddArray(arrayName);
    if(allContents[idx] == 0)
      {
      s->DisableArray(arrayName);
      }
    else
      {
      s->EnableArray(arrayName);
      }
    }
  delete [] allContents;
}

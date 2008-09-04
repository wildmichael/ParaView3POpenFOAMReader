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

#include "pqPOpenFOAMPanel.h"

// Qt
#include <QComboBox>
#include <QLayout>

// server manager
#include "vtkSMIntVectorProperty.h"
#include "vtkSMSourceProxy.h"

//-----------------------------------------------------------------------------
pqPOpenFOAMPanel::pqPOpenFOAMPanel(pqProxy *pxy, QWidget *p)
  : pqOpenFOAMPanel(pxy, p)
{
  // create case selection combo box
  QComboBox *type = new QComboBox();
  type->addItem("Decomposed case");
  type->addItem("Reconstructed case");
  // place the box at top-left of the layout grid
  this->PanelLayout->addWidget(type, 0, 0);
  vtkSMIntVectorProperty::SafeDownCast(vtkSMSourceProxy::SafeDownCast(
    this->proxy())->GetProperty("CaseType"))->SetImmediateUpdate(1);
  QObject::connect(type, SIGNAL(currentIndexChanged(int)), this,
    SLOT(onCurrentIndexChanged(int)));
}

//-----------------------------------------------------------------------------
void pqPOpenFOAMPanel::onCurrentIndexChanged(const int index)
{
  vtkSMSourceProxy* sp = vtkSMSourceProxy::SafeDownCast(this->proxy());
  vtkSMIntVectorProperty *propCaseType
    = vtkSMIntVectorProperty::SafeDownCast(sp->GetProperty("CaseType"));
  propCaseType->SetElements1(index);
  // force updating everything (RequestInformation() + RequestData())
  sp->UpdatePipeline();
}

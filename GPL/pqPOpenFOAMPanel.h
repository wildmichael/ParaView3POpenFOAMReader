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

#ifndef _pqPOpenFOAMPanel_h
#define _pqPOpenFOAMPanel_h

#include "pqOpenFOAMPanel.h"

class pqPOpenFOAMPanel : public pqOpenFOAMPanel
{
  Q_OBJECT
public:
  pqPOpenFOAMPanel(pqProxy *, QWidget *);

  ~pqPOpenFOAMPanel() {}

protected slots:
  void onCurrentIndexChanged(const int);
};

#endif
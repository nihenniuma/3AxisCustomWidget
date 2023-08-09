#include "vtkAxisWidget.h"

vtkStandardNewMacro(vtkAxisWidget);

void vtkAxisWidget::PlaceWidget(double bds[6]) {
  double bounds[6];
  double center[3];
  this->AdjustBounds(bds,bounds,center);
  // this->PositionHandles();
  this->SizeHandles();
}

vtkAxisWidget::vtkAxisWidget() {
  this->_XCircle = vtkRegularPolygonSource::New();
  this->_XCircleMapper = vtkPolyDataMapper::New();
  this->_XCircleActor = vtkActor::New();
  this->_XCircle->SetNumberOfSides(100);
  this->_XCircle->SetCenter(0, 0, 0);
  this->_XCircle->SetRadius(100);
  this->_XCircleMapper->SetInputData(_XCircle->GetOutput());
  this->_XCircleActor->SetMapper(this->_XCircleMapper);
  this->_XCircleActor->GetProperty()->SetColor(1, 1, 0);
}

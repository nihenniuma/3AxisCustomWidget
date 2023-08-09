#ifndef VTK_AXIS_WIDGET_H
#define VTK_AXIS_WIDGET_H

#include <vtkWidgetRepresentation.h>
#include <vtkAbstractWidget.h>
#include <vtkWidgetRepresentation.h>
#include <vtkActor.h>
#include <vtkRegularPolygonSource.h>
#include <vtkMapper.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkConeSource.h>
#include <vtkCellPicker.h>
#include <vtkLineSource.h>
#include <vtkAssembly.h>
#include <vtkMatrix4x4.h>
#include <vtkTransform.h>
#include <vtkPoints.h>
#include <vtkDoubleArray.h>
#include <vtkRenderer.h>
#include <vtkCamera.h>
#include <vtkWindow.h>
#include <vtkInteractionWidgetsModule.h>
#include <vtk3DWidget.h>

class VTKINTERACTIONWIDGETS_EXPORT vtkAxisWidget : vtk3DWidget{
public:
  static vtkAxisWidget* New();

  vtkTypeMacro(vtkAxisWidget,vtk3DWidget);

  void PlaceWidget (double bounds[6]);


protected:
  vtkAxisWidget();
  // ~vtkAxisWidget();

private:
  vtkRegularPolygonSource* _XCircle = nullptr;
  vtkPolyDataMapper* _XCircleMapper = nullptr;
  vtkActor* _XCircleActor = nullptr;
};
#endif
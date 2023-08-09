#ifndef VTK_CUSTOM_REPRESENTATION_H
#define VTK_CUSTOM_REPRESENTATION_H

#include <vtkWidgetRepresentation.h>
#include <vtkAbstractWidget.h>
#include <vtkWidgetRepresentation.h>
#include <vtkActor.h>
#include <vtkRegularPolygonSource.h>
#include <vtkLineSource.h>
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
#include <vtkSphereSource.h>

class vtkArrowSource;
class vtkConeSource;
class vtkAppendPolyData;

class vtkCustomRepresentation : public vtkWidgetRepresentation{
public:
  static vtkCustomRepresentation* New();

  vtkTypeMacro(vtkCustomRepresentation, vtkWidgetRepresentation);

  void BuildRepresentation() override;
  void SizeHandles();
  void SetInteractionState(int state);
  int ComputeInteractionState(int X, int Y, int modify);
  void StartWidgetInteraction(double e[2]) override;
  void WidgetInteraction(double e[2]) override;
  void StartComplexInteraction(
    vtkRenderWindowInteractor *iren,
    vtkAbstractWidget *widget,
    unsigned long event, void *calldata) override;
  void ComplexInteraction(
    vtkRenderWindowInteractor *iren,
    vtkAbstractWidget *widget,
    unsigned long event, void *calldata) override;
  // int ComputeComplexInteractionState(
  //   vtkRenderWindowInteractor *iren,
  //   vtkAbstractWidget *widget,
  //   unsigned long event, void *calldata, int modify = 0) override;
  // void EndComplexInteraction(
  //   vtkRenderWindowInteractor *iren,
  //   vtkAbstractWidget *widget,
  //   unsigned long event, void *calldata) override;

  // double* GetBounds() VTK_SIZEHINT(6) override;
  void PlaceWidget(double bounds[6]) override;
  // void GetActors(vtkPropCollection* pc) override;
  void ReleaseGraphicsResources(vtkWindow*) override;
  int RenderOpaqueGeometry(vtkViewport*) override;
  int RenderTranslucentPolygonalGeometry(vtkViewport*) override;
  // vtkTypeBool HasTranslucentPolygonalGeometry() override;

  void Rotate(double *p1, double *p2, int singleAxis);
  void RotateByMatrix(vtkMatrix4x4 *origin_matrix, double *p1, double *p2,int direction,vtkMatrix4x4* result_matrix);
  void GetPlaneProjection(double*, double*, double*);
  void MoveAxis(const int singleAxis, double *p1, double *p2);
  void MoveByMatrix(vtkMatrix4x4 *origin_matrix, double *p1, double *p2,int direction,vtkMatrix4x4* result_matrix);

  void SetTargetActor(vtkActor** actor);
  virtual void CreateDefaultProperties(); 

  enum {Outside = 0, RotateAroundX, RotateAroundY, RotateAroundZ, MoveAlongX, MoveAlongY, MoveAlongZ, Moving};
protected:
  vtkCustomRepresentation();
  // ~vtkCustomRepresentation();

  int Highlight(vtkProp* prop);
  void HightlightCircle(int);
  void HightlightAxis(int);

  vtkActor* CurrentActor = nullptr;
  double LastEventPosition[3];
  double LastEventOrientation[4];
  double StartEventOrientation[4];
  double SnappedEventOrientations[3][4];
  bool SnappedOrientation[3];
  bool SnapToAxes;
private:
  // target actor
  vtkActor* TargetActor = nullptr;

  // Axis Line
  vtkLineSource** AxisLineSource = nullptr;
  vtkConeSource** AxisConeSource = nullptr;
  vtkAppendPolyData** AxisWithArrow = nullptr;
  vtkPolyDataMapper** AxisWithArrowMapper = nullptr;
  vtkActor** Axis = nullptr;

  vtkActor* CenterActor = nullptr;
  vtkPolyDataMapper* CenterMapper = nullptr;
  vtkSphereSource* CenterSource = nullptr;

  vtkActor** AxisCircle = nullptr;
  vtkPolyDataMapper** AxisCircleMapper = nullptr;
  vtkRegularPolygonSource** AxisPolygonSource = nullptr;
  double AxisColor[3][3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
  double AxisLineWidth = 2;
  double AxisLineLength = 2;

  vtkProperty** CircleProperty = nullptr;
  vtkProperty** SelectedCircleProperty = nullptr;
  vtkProperty** AxisProperty = nullptr;
  vtkProperty** SelectedAxisProperty = nullptr;
  vtkProperty* CenterProperty = nullptr;
  vtkProperty* SelectedCenterProperty = nullptr;

  vtkCellPicker* Picker = nullptr;
  vtkCellPicker* LastPicker = nullptr;

  vtkMatrix4x4   *Matrix;
  // vtkCellPicker* CirclePicker = nullptr;
  // vtkCellPicker* AxisPicker = nullptr;

  /*-------------test1---------------------*/
  // vtkCellPicker* CirclePicker = nullptr;
  // vtkCellPicker* AxisPicker = nullptr;
  // vtkTransform* Transform = nullptr;
  // vtkPoints* PlanePoints = nullptr;
  // vtkDoubleArray* PlaneNormals = nullptr;
  // vtkMatrix4x4* Matrix = nullptr;

  // vtkActor* AxisCircle[3];
  // vtkActor* AxisLine[3];
  // vtkRegularPolygonSource* AxisPolygonSource[3];
  // vtkPolyDataMapper* AxisCircleMapper[3];
  // vtkAssembly* HandleAxis[3];

  // double g_normalColor[3][3] = {1.0};
  // double g_circleRadius = 10;
  // double g_circleLineWidth = 1;
  // double g_axisLength = 5;
  // double g_axisLineWidth = 1;
  /*-------------test1---------------------*/
  // vtkRegularPolygonSource* _XCircle = nullptr;
  // vtkRegularPolygonSource* _YCircle = nullptr;
  // vtkRegularPolygonSource* _ZCircle = nullptr;

  // vtkPolyDataMapper* _XCircleMapper = nullptr;
  // vtkActor* _XCircleActor = nullptr;
  /*-------------test2---------------------*/
};

#endif
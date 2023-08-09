#ifndef VTK_CUSTOM_WIDGET_H
#define VTK_CUSTOM_WIDGET_H

#include "vtkCustomRepresentation.h"

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
#include <vtkInteractionWidgetsModule.h>

class vtkWidgetCallbackMapper;

class vtkCustomWidget : public vtkAbstractWidget{
public:
  static vtkCustomWidget* New();

  vtkTypeMacro(vtkCustomWidget, vtkAbstractWidget);

  void SetRepresentation(vtkCustomRepresentation* r)
  {
    this->Superclass::SetWidgetRepresentation(reinterpret_cast<vtkWidgetRepresentation*>(r));
  }

  void SetTargetActor(vtkActor** actor);

  void CreateDefaultRepresentation();
protected:
  vtkCustomWidget();
  // ~vtkCustomWidget();

  vtkSetMacro(TranslationEnabled,vtkTypeBool);
  vtkGetMacro(TranslationEnabled,vtkTypeBool);
  vtkBooleanMacro(TranslationEnabled,vtkTypeBool);
  vtkSetMacro(ScalingEnabled,vtkTypeBool);
  vtkGetMacro(ScalingEnabled,vtkTypeBool);
  vtkBooleanMacro(ScalingEnabled,vtkTypeBool);
  vtkSetMacro(RotationEnabled,vtkTypeBool);
  vtkGetMacro(RotationEnabled,vtkTypeBool);
  vtkBooleanMacro(RotationEnabled,vtkTypeBool);
  vtkSetMacro(CircleRotationEnabled,vtkTypeBool);
  vtkGetMacro(CircleRotationEnabled,vtkTypeBool);
  vtkBooleanMacro(CircleRotationEnabled,vtkTypeBool);
  vtkSetMacro(LineTranslationEnabled,vtkTypeBool);
  vtkGetMacro(LineTranslationEnabled,vtkTypeBool);
  vtkBooleanMacro(LineTranslationEnabled,vtkTypeBool);

  static void SelectAction(vtkAbstractWidget*);
  static void EndSelectAction(vtkAbstractWidget*);
  static void RotateAction(vtkAbstractWidget*);
  static void MoveAxisAction(vtkAbstractWidget*);
  static void RotateAction3D(vtkAbstractWidget*);
  static void SelectAction3D(vtkAbstractWidget *w);
  static void EndSelectAction3D(vtkAbstractWidget *w);

  int WidgetState = 0;
  enum _WidgetState {Start = 0, Active};

  vtkTypeBool TranslationEnabled;
  vtkTypeBool ScalingEnabled;
  vtkTypeBool RotationEnabled;
  vtkTypeBool CircleRotationEnabled;
  vtkTypeBool LineTranslationEnabled;

private:

};
#endif
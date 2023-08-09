#include "vtkCustomRepresentation.h"
#include <iostream>
#include <vtkAssemblyPath.h>
#include <vtkArrowSource.h>
#include <vtkConeSource.h>
#include <vtkEventData.h>
#include <vtkAppendPolyData.h>

vtkStandardNewMacro(vtkCustomRepresentation);

vtkCustomRepresentation::vtkCustomRepresentation() {
  this->InteractionState = vtkCustomRepresentation::Outside;

  this->CreateDefaultProperties();

  // Center
  this->CenterActor = vtkActor::New();
  this->CenterMapper = vtkPolyDataMapper::New();
  this->CenterSource = vtkSphereSource::New();
  this->CenterSource->SetThetaResolution(16);
  this->CenterSource->SetPhiResolution(16);
  this->CenterSource->SetRadius(0.25);
  this->CenterMapper->SetInputConnection(this->CenterSource->GetOutputPort());
  this->CenterActor->SetMapper(this->CenterMapper);
  this->CenterActor->SetProperty(this->CenterProperty);
  vtkNew<vtkMatrix4x4> initMatrix;
  this->CenterActor->SetUserMatrix(initMatrix); // if we want to use GetUserMatrix, this step is needed!

  // Circle
  this->AxisCircle = new vtkActor* [3];
  this->AxisPolygonSource = new vtkRegularPolygonSource* [3];
  this->AxisCircleMapper = new vtkPolyDataMapper* [3];
  for(int i = 0; i < 3; ++i){
    // Circle Axis
    this->AxisPolygonSource[i] = vtkRegularPolygonSource::New();
    this->AxisPolygonSource[i]->SetNumberOfSides(360);
    this->AxisPolygonSource[i]->SetCenter(0, 0, 0);
    this->AxisPolygonSource[i]->SetRadius(1);

    this->AxisCircleMapper[i] = vtkPolyDataMapper::New();
    this->AxisCircleMapper[i]->SetInputConnection(this->AxisPolygonSource[i]->GetOutputPort());
    
    this->AxisCircle[i] = vtkActor::New();
    this->AxisCircle[i]->SetMapper(this->AxisCircleMapper[i]);
    this->AxisCircle[i]->SetProperty(this->CircleProperty[i]);
    if(i == 0){
      this->AxisPolygonSource[i]->SetNormal(1,0,0);
    }else if(i == 1){
      this->AxisPolygonSource[i]->SetNormal(0,1,0);
    }else if(i == 2){
      this->AxisPolygonSource[i]->SetNormal(0,0,1);
    }

    vtkNew<vtkMatrix4x4> initMatrix;
    this->AxisCircle[i]->SetUserMatrix(initMatrix);
  }

  // Axis
  this->Axis = new vtkActor* [3];
  this->AxisWithArrow = new vtkAppendPolyData* [3];
  this->AxisWithArrowMapper = new vtkPolyDataMapper* [3];
  this->AxisConeSource = new vtkConeSource* [3];
  this->AxisLineSource = new vtkLineSource* [3];
  for(int i = 0; i < 3; ++i){
    double Point1[3] = {0};
    Point1[i] = this->AxisLineLength;
    double Point2[3] = {0};
    Point2[i] = -this->AxisLineLength;
    this->AxisLineSource[i] = vtkLineSource::New();
    this->AxisLineSource[i]->SetPoint1(Point1);
    this->AxisLineSource[i]->SetPoint2(Point2);

    double Dir[3] = {0};
    Dir[i] = this->AxisLineLength;
    this->AxisConeSource[i] = vtkConeSource::New();
    this->AxisConeSource[i]->SetHeight(0.3);
    this->AxisConeSource[i]->SetRadius(0.1);
    this->AxisConeSource[i]->SetResolution(20);
    this->AxisConeSource[i]->SetDirection(Dir);
    this->AxisConeSource[i]->SetCenter(Dir);

    this->AxisWithArrow[i] = vtkAppendPolyData::New();
    this->AxisWithArrow[i]->AddInputConnection(this->AxisConeSource[i]->GetOutputPort());
    this->AxisWithArrow[i]->AddInputConnection(this->AxisLineSource[i]->GetOutputPort());
    
    this->AxisWithArrowMapper[i] = vtkPolyDataMapper::New();
    this->AxisWithArrowMapper[i]->SetInputConnection(this->AxisWithArrow[i]->GetOutputPort());

    this->Axis[i] = vtkActor::New();
    this->Axis[i]->SetMapper(AxisWithArrowMapper[i]);
    this->Axis[i]->SetProperty(this->AxisProperty[i]);
    this->Axis[i]->PickableOn();

    vtkNew<vtkMatrix4x4> initMatrix;
    this->Axis[i]->SetUserMatrix(initMatrix);
  }

  this->Picker = vtkCellPicker::New();
  for(int i = 0; i < 3; ++i){
    this->Picker->AddPickList(this->AxisCircle[i]);
  }
  for(int i = 0; i < 3; ++i){
    this->Picker->AddPickList(this->Axis[i]);
  }
  this->Picker->AddPickList(this->CenterActor);
  this->Picker->SetTolerance(0.001);
  this->Picker->PickFromListOn();

  this->Matrix = vtkMatrix4x4::New();
}

int vtkCustomRepresentation::Highlight(vtkProp* prop) { 
  if(this->CurrentActor){
    if(this->CurrentActor == this->CenterActor){
      std::cout << "delight center" << std::endl;
      this->CurrentActor->SetProperty(this->CenterProperty);
    }else{
      for(int i = 0; i < 3; ++i){
        if(this->AxisCircle[i] == this->CurrentActor){
          std::cout << "delight rotate" << std::endl;
          this->CurrentActor->SetProperty(this->CircleProperty[i]);
          break;
        }else if(this->Axis[i] == this->CurrentActor){
          std::cout << "delight axis" << std::endl;
          this->Axis[i]->SetProperty(this->AxisProperty[i]);
          break;
        }else{
          continue;
        }
      }
    }
  }
  this->CurrentActor = static_cast<vtkActor *>(prop);
  if(this->CurrentActor){
    if(this->CurrentActor == this->CenterActor){
      std::cout << "highlight center" << std::endl;
      this->CurrentActor->SetProperty(this->SelectedCenterProperty);
    }else{
      for(int i = 0; i < 3; ++i){
        if(this->AxisCircle[i] == this->CurrentActor){
          std::cout << "highlight rotate" << std::endl;
          this->CurrentActor->SetProperty(this->SelectedCircleProperty[i]);
          return i;
        }else if(this->Axis[i] == this->CurrentActor){
          std::cout << "highlight axis" << std::endl;
          this->Axis[i]->SetProperty(this->SelectedAxisProperty[i]);
          return i;
        }else{
          continue;
        }
      }
    }
  }
  return -1;
}

void vtkCustomRepresentation::HightlightCircle(int index) {}

void vtkCustomRepresentation::HightlightAxis(int index) {}

void vtkCustomRepresentation::BuildRepresentation() {
  if (this->GetMTime() > this->BuildTime ||
    (this->Renderer && this->Renderer->GetVTKWindow() &&
      (this->Renderer->GetVTKWindow()->GetMTime() > this->BuildTime ||
        this->Renderer->GetActiveCamera()->GetMTime() > this->BuildTime)))
  {
    this->SizeHandles();
    this->BuildTime.Modified();
  }
}

void vtkCustomRepresentation::SizeHandles() {

}

void vtkCustomRepresentation::SetInteractionState(int state) {
  state = (state < vtkCustomRepresentation::Outside ? vtkCustomRepresentation::Outside : 
            (state > vtkCustomRepresentation::Moving ? vtkCustomRepresentation::Moving : state));

  this->InteractionState = state;
  switch (state){
    case vtkCustomRepresentation::RotateAroundX:
      Highlight(this->AxisCircle[0]);
      break;
    case vtkCustomRepresentation::RotateAroundY:
      Highlight(this->AxisCircle[1]);
      break;
    case vtkCustomRepresentation::RotateAroundZ:
      Highlight(this->AxisCircle[2]);
      break;
    case vtkCustomRepresentation::MoveAlongX:
      Highlight(this->Axis[0]);
      break;
    case vtkCustomRepresentation::MoveAlongY:
      Highlight(this->Axis[1]);
      break;
    case vtkCustomRepresentation::MoveAlongZ:
      Highlight(this->Axis[2]);
      break;
    case vtkCustomRepresentation::Moving:
      Highlight(this->CenterActor);
      break;
    default:
      Highlight(nullptr);
      break;

  }
}

int vtkCustomRepresentation::ComputeInteractionState(int X, int Y, int modify) {
  if (!this->Renderer || !this->Renderer->IsInViewport(X, Y))
  {
    this->InteractionState = vtkCustomRepresentation::Outside;
    return this->InteractionState;
  }
  
  vtkAssemblyPath* path = this->GetAssemblyPath(X, Y, 0., this->Picker);
  if ( path != nullptr )
  {
    this->ValidPick = 1;
    this->LastPicker = this->Picker;
    this->CurrentActor =
           reinterpret_cast<vtkActor *>(path->GetFirstNode()->GetViewProp());
    if ( this->CurrentActor == this->AxisCircle[0] )
    {
      std::cout << "RotateAroundX" << std::endl;
      this->InteractionState = vtkCustomRepresentation::RotateAroundX;
    }
    else if ( this->CurrentActor == this->AxisCircle[1] )
    {
      std::cout << "RotateAroundY" << std::endl;
      this->InteractionState = vtkCustomRepresentation::RotateAroundY;
    }
    else if ( this->CurrentActor == this->AxisCircle[2] )
    {
      std::cout << "RotateAroundZ" << std::endl;
      this->InteractionState = vtkCustomRepresentation::RotateAroundZ;
    }
    else if ( this->CurrentActor == this->Axis[0] )
    {
      std::cout << "MoveAlongX" << std::endl;
      this->InteractionState = vtkCustomRepresentation::MoveAlongX;
    }
    else if ( this->CurrentActor == this->Axis[1] )
    {
      std::cout << "MoveAlongY" << std::endl;
      this->InteractionState = vtkCustomRepresentation::MoveAlongY;
    }
    else if ( this->CurrentActor == this->Axis[2] )
    {
      std::cout << "MoveAlongZ" << std::endl;
      this->InteractionState = vtkCustomRepresentation::MoveAlongZ;
    }
    else if ( this->CurrentActor == this->CenterActor )
    {
      std::cout << "Move" << std::endl;
      this->InteractionState = vtkCustomRepresentation::Moving;
    }
  }else{
    this->InteractionState = vtkCustomRepresentation::Outside;
  }

  return this->InteractionState;
}

void vtkCustomRepresentation::StartWidgetInteraction(double e[2]) {
  // Store the start position
  this->StartEventPosition[0] = e[0];
  this->StartEventPosition[1] = e[1];
  this->StartEventPosition[2] = 0.0;

  // Store the start position
  this->LastEventPosition[0] = e[0];
  this->LastEventPosition[1] = e[1];
  this->LastEventPosition[2] = 0.0;

  this->ComputeInteractionState(static_cast<int>(e[0]),static_cast<int>(e[1]),0);
}

void vtkCustomRepresentation::WidgetInteraction(double e[2]) {
  // Convert events to appropriate coordinate systems
  vtkCamera *camera = this->Renderer->GetActiveCamera();
  if ( !camera )
  {
    return;
  }
  double focalPoint[4], pickPoint[4], prevPickPoint[4];
  double z, vpn[3];
  camera->GetViewPlaneNormal(vpn);

  // Compute the two points defining the motion vector
  double pos[3];
  this->Picker->GetPickPosition(pos);

  vtkInteractorObserver::ComputeWorldToDisplay(this->Renderer,
                                               pos[0], pos[1], pos[2],
                                               focalPoint);
  z = focalPoint[2];
  vtkInteractorObserver::ComputeDisplayToWorld(this->Renderer,this->LastEventPosition[0],
                                               this->LastEventPosition[1], z, prevPickPoint);
  vtkInteractorObserver::ComputeDisplayToWorld(this->Renderer, e[0], e[1], z, pickPoint);

  // Process the motion
  if ( this->InteractionState == vtkCustomRepresentation::RotateAroundX )
  {
    std::cout << "start rotating around z!" << std::endl;
    this->Rotate(prevPickPoint, pickPoint, 0);
  }

  else if ( this->InteractionState == vtkCustomRepresentation::RotateAroundY )
  {
    std::cout << "start rotating around z!" << std::endl;
    this->Rotate(prevPickPoint, pickPoint, 1);
  }

  else if ( this->InteractionState == vtkCustomRepresentation::RotateAroundZ )
  {
    std::cout << "start rotating around z!" << std::endl;
    this->Rotate(prevPickPoint, pickPoint, 2);
  }

  else if ( this->InteractionState == vtkCustomRepresentation::MoveAlongX )
  {
    std::cout << "start moving along x!" << std::endl;
    this->MoveAxis(0, prevPickPoint, pickPoint);
  }

  else if ( this->InteractionState == vtkCustomRepresentation::MoveAlongY )
  {
    std::cout << "start moving along y!" << std::endl;
    this->MoveAxis(1, prevPickPoint, pickPoint);
  }

  else if ( this->InteractionState == vtkCustomRepresentation::MoveAlongZ )
  {
    std::cout << "start moving along z!" << std::endl;
    this->MoveAxis(2, prevPickPoint, pickPoint);
  }

  else if ( this->InteractionState == vtkCustomRepresentation::Moving )
  {
    std::cout << "start moving!" << std::endl;
    // this->UpdatePose(this->LastEventPosition, this->LastEventOrientation,
    //                   eventPos, eventDir);
  }

  // Store the start position
  this->LastEventPosition[0] = e[0];
  this->LastEventPosition[1] = e[1];
  this->LastEventPosition[2] = 0.0;
}

void vtkCustomRepresentation::StartComplexInteraction(
  vtkRenderWindowInteractor* iren, vtkAbstractWidget* widget,
  unsigned long event, void* calldata) 
{
  vtkEventData *edata = static_cast<vtkEventData *>(calldata);
  vtkEventDataDevice3D *edd = edata->GetAsEventDataDevice3D();
  if (edd)
  {
    edd->GetWorldPosition(this->StartEventPosition);
    this->LastEventPosition[0] = this->StartEventPosition[0];
    this->LastEventPosition[1] = this->StartEventPosition[1];
    this->LastEventPosition[2] = this->StartEventPosition[2];
    edd->GetWorldOrientation(this->StartEventOrientation);
    std::copy(this->StartEventOrientation, this->StartEventOrientation + 4,
      this->LastEventOrientation);
    for (int i = 0; i < 3; ++i)
    {
      if (this->SnappedOrientation[i])
      {
        std::copy(this->StartEventOrientation,
          this->StartEventOrientation+4, this->SnappedEventOrientations[i]);
      }
    }
  }
}

void vtkCustomRepresentation::ComplexInteraction(vtkRenderWindowInteractor*,
                                                 vtkAbstractWidget*,
                                                 unsigned long,
                                                 void* calldata) {
  vtkEventData *edata = static_cast<vtkEventData *>(calldata);
  vtkEventDataDevice3D *edd = edata->GetAsEventDataDevice3D();
  if (edd)
  {
    // all others
    double eventPos[3];
    edd->GetWorldPosition(eventPos);
    double eventDir[4];
    edd->GetWorldOrientation(eventDir);

    double *prevPickPoint = this->LastEventPosition;
    double *pickPoint = eventPos;

    if ( this->InteractionState == vtkCustomRepresentation::RotateAroundX )
    {
      std::cout << "start rotating around z!" << std::endl;
      this->Rotate(prevPickPoint, pickPoint, 0);
    }

    else if ( this->InteractionState == vtkCustomRepresentation::RotateAroundY )
    {
      std::cout << "start rotating around z!" << std::endl;
      this->Rotate(prevPickPoint, pickPoint, 1);
    }

    else if ( this->InteractionState == vtkCustomRepresentation::RotateAroundZ )
    {
      std::cout << "start rotating around z!" << std::endl;
      this->Rotate(prevPickPoint, pickPoint, 2);
    }

    else if ( this->InteractionState == vtkCustomRepresentation::MoveAlongX )
    {
      std::cout << "start moving along x!" << std::endl;
      // this->MovePlusYFace(prevPickPoint,pickPoint);
    }

    else if ( this->InteractionState == vtkCustomRepresentation::MoveAlongY )
    {
      std::cout << "start moving along y!" << std::endl;
      // this->MoveMinusZFace(prevPickPoint,pickPoint);
    }

    else if ( this->InteractionState == vtkCustomRepresentation::MoveAlongZ )
    {
      std::cout << "start moving along z!" << std::endl;
      // this->MovePlusZFace(prevPickPoint,pickPoint);
    }

    else if ( this->InteractionState == vtkCustomRepresentation::Moving )
    {
      std::cout << "start moving!" << std::endl;
      // this->UpdatePose(this->LastEventPosition, this->LastEventOrientation,
      //                   eventPos, eventDir);
    }

    // Book keeping
    std::copy(eventPos, eventPos + 3, this->LastEventPosition);
    std::copy(eventDir, eventDir + 4, this->LastEventOrientation);
    this->Modified();
  }                                              
}

void vtkCustomRepresentation::PlaceWidget(double bds[6]) {
  int i;
  double bounds[6], center[3];

  this->AdjustBounds(bds,bounds,center);

  // this->Points->SetPoint(0, bounds[0], bounds[2], bounds[4]);
  // this->Points->SetPoint(1, bounds[1], bounds[2], bounds[4]);
  // this->Points->SetPoint(2, bounds[1], bounds[3], bounds[4]);
  // this->Points->SetPoint(3, bounds[0], bounds[3], bounds[4]);
  // this->Points->SetPoint(4, bounds[0], bounds[2], bounds[5]);
  // this->Points->SetPoint(5, bounds[1], bounds[2], bounds[5]);
  // this->Points->SetPoint(6, bounds[1], bounds[3], bounds[5]);
  // this->Points->SetPoint(7, bounds[0], bounds[3], bounds[5]);

  // for (i=0; i<6; i++)
  // {
  //   this->InitialBounds[i] = bounds[i];
  // }
  // this->InitialLength = sqrt((bounds[1]-bounds[0])*(bounds[1]-bounds[0]) +
  //                            (bounds[3]-bounds[2])*(bounds[3]-bounds[2]) +
  //                            (bounds[5]-bounds[4])*(bounds[5]-bounds[4]));

  // this->PositionHandles();
  // this->ComputeNormals();
  this->ValidPick = 1; //since we have set up widget
  this->SizeHandles();
}

void vtkCustomRepresentation::ReleaseGraphicsResources(vtkWindow*) {}

int vtkCustomRepresentation::RenderOpaqueGeometry(vtkViewport* v) {
  int count = 0;
  this->BuildRepresentation();

  count += this->CenterActor->RenderOpaqueGeometry(v);
  for(int i = 0; i < 3; i++){
    count += this->AxisCircle[i]->RenderOpaqueGeometry(v);
    // count += this->HandleAxis[i]->RenderOpaqueGeometry(v);
    count += this->Axis[i]->RenderOpaqueGeometry(v);
  }
  if(this->TargetActor){
    std::cout << "render target! \n";
    count += this->TargetActor->RenderOpaqueGeometry(v);
  }
  return count;
}

int vtkCustomRepresentation::RenderTranslucentPolygonalGeometry(vtkViewport* v) {
  int count = 0;
  this->BuildRepresentation();

  count += this->CenterActor->RenderTranslucentPolygonalGeometry(v);
  for(int i = 0; i < 3; i++){
    count += this->AxisCircle[i]->RenderTranslucentPolygonalGeometry(v);
    // count += this->HandleAxis[i]->RenderTranslucentPolygonalGeometry(v);
    count += this->Axis[i]->RenderTranslucentPolygonalGeometry(v);
  }
  if(this->TargetActor){
    std::cout << "render target! \n";
    count += this->TargetActor->RenderTranslucentPolygonalGeometry(v);
  }
  return count;
}

void vtkCustomRepresentation::Rotate(double *p1, double *p2, int singleAxis)
{
    for(int i=0;i<3;i++){
      vtkMatrix4x4 *origin_matrixCircle = this->AxisCircle[i]->GetUserMatrix();
      vtkNew<vtkMatrix4x4> result_matrixCircle;
      RotateByMatrix(origin_matrixCircle,p1,p2,singleAxis,result_matrixCircle);
      std::cout << "AxisCircle Matrix calculate!" << std::endl;
      for(int i = 0; i < 4; ++i){
        for(int j = 0; j < 4; ++j){
          std::cout << result_matrixCircle->GetElement(i, j) << " ";
        }
        std::cout << std::endl;
      }
      this->AxisCircle[i]->SetUserMatrix(result_matrixCircle);

      vtkMatrix4x4 *origin_matrixAxis = this->Axis[i]->GetUserMatrix();
      vtkNew<vtkMatrix4x4> result_matrixAxis;
      RotateByMatrix(origin_matrixAxis,p1,p2,singleAxis,result_matrixAxis);
      std::cout << "Axis Matrix calculate!" << std::endl;
      this->Axis[i]->SetUserMatrix(result_matrixAxis);
    }
    if(this->TargetActor){
      vtkMatrix4x4 *origin_matrix = this->TargetActor->GetUserMatrix();
      vtkNew<vtkMatrix4x4> result_matrix;
      RotateByMatrix(origin_matrix,p1,p2,singleAxis,result_matrix);
      std::cout << "TargetActor Matrix calculate!" << std::endl;
      for(int i = 0; i < 4; ++i){
        for(int j = 0; j < 4; ++j){
          std::cout << result_matrix->GetElement(i, j) << " ";
        }
        std::cout << std::endl;
      }
      this->TargetActor->SetUserMatrix(result_matrix);
    }

    this->Modified();
    this->BuildRepresentation();
}

void vtkCustomRepresentation::RotateByMatrix(vtkMatrix4x4 *origin_matrix, double *p1, double *p2,int direction,vtkMatrix4x4* result_matrix)
{
  // 将鼠标位置移动到自身坐标系下，求两次鼠标位置向量在投影平面的夹角
  vtkNew<vtkTransform> trans;
  trans->SetMatrix(origin_matrix);
  double pos_t1[4] { p1[0], p1[1], p1[2], 1 };
  double pos_t2[4] { p2[0], p2[1], p2[2], 1 };
  vtkNew<vtkMatrix4x4> posture_inv;
  vtkMatrix4x4::Invert(origin_matrix, posture_inv);
  auto pos_t = posture_inv->MultiplyDoublePoint(pos_t1);
  double v1[3] = { pos_t[0], pos_t[1], pos_t[2] };
  pos_t = posture_inv->MultiplyDoublePoint(pos_t2);
  double v2[3] = { pos_t[0], pos_t[1], pos_t[2] };
  double normal[3];

  if(direction==0){
    normal[0] = 1;
    normal[1] = 0;
    normal[2] = 0;
  }else if(direction==1){
    normal[0] = 0;
    normal[1] = 1;
    normal[2] = 0;
  }else if(direction==2){
    normal[0] = 0;
    normal[1] = 0;
    normal[2] = 1;
  }else{
    return;
  }
  double projection1[3], projection2[3];

  GetPlaneProjection(v1, normal, projection1);
  std::cout << "projection1: \n";
  for(int i = 0 ; i < 3; ++i){
    std::cout << projection1[i] << " ";
  }
  std::cout << "\n";
  GetPlaneProjection(v2, normal, projection2);
  std::cout << "projection2: \n";
  for(int i = 0 ; i < 3; ++i){
    std::cout << projection2[i] << " ";
  }
  std::cout << "\n";

  vtkMath::Normalize(projection1);
  vtkMath::Normalize(projection2);
  double axis[3];
  vtkMath::Cross(projection1, projection2, axis);
  double radians = acos(vtkMath::Dot(projection1, projection2));
  double degrees = vtkMath::DegreesFromRadians(radians);
  trans->RotateWXYZ(degrees, axis);
  result_matrix->DeepCopy(trans->GetMatrix());
}

void vtkCustomRepresentation::MoveAxis(const int singleAxis, double *p1, double *p2)
{
  // for(int i = 0; i < 3; ++i){
  //   vtkMatrix4x4 *origin_matrixCircle = this->AxisCircle[i]->GetUserMatrix();
  //   vtkNew<vtkMatrix4x4> result_matrixCircle;
  //   MoveByMatrix(origin_matrixCircle,p1,p2,singleAxis,result_matrixCircle);
  //   std::cout << "AxisCircle Matrix calculate!" << std::endl;
  //   for(int i = 0; i < 4; ++i){
  //     for(int j = 0; j < 4; ++j){
  //       std::cout << result_matrixCircle->GetElement(i, j) << " ";
  //     }
  //     std::cout << std::endl;
  //   }
  //   this->AxisCircle[i]->SetUserMatrix(result_matrixCircle);

  //   vtkMatrix4x4* origin_matrixAxis = this->Axis[i]->GetUserMatrix();
  //   vtkNew<vtkMatrix4x4> result_matrixAxis;
  //   MoveByMatrix(origin_matrixAxis,p1,p2,singleAxis,result_matrixAxis);
  //   std::cout << "AxisCircle Matrix calculate!" << std::endl;
  //   for(int i = 0; i < 4; ++i){
  //     for(int j = 0; j < 4; ++j){
  //       std::cout << result_matrixAxis->GetElement(i, j) << " ";
  //     }
  //     std::cout << std::endl;
  //   }
  //   this->Axis[i]->SetUserMatrix(result_matrixAxis);
  // }
  
  vtkNew<vtkMatrix4x4> result_matrix;
  for(int i = 0; i < 3; ++i){
    vtkMatrix4x4* origin_matrix = this->Axis[i]->GetUserMatrix();
    MoveByMatrix(origin_matrix,p1,p2,singleAxis,result_matrix);
    std::cout << "AxisCircle Matrix calculate!" << std::endl;
    for(int i = 0; i < 4; ++i){
      for(int j = 0; j < 4; ++j){
        std::cout << result_matrix->GetElement(i, j) << " ";
      }
      std::cout << std::endl;
    }
    this->Axis[i]->SetUserMatrix(result_matrix);
    this->AxisCircle[i]->SetUserMatrix(result_matrix);
  }
  this->CenterActor->SetUserMatrix(result_matrix);
  if(this->TargetActor){
    this->TargetActor->SetUserMatrix(result_matrix);
  }

  this->Modified();
  this->BuildRepresentation();
}

void vtkCustomRepresentation::MoveByMatrix(vtkMatrix4x4* origin_matrix,
                                           double* p1, double* p2,
                                           int direction,
                                           vtkMatrix4x4* result_matrix) 
{
  vtkNew<vtkTransform> posture_transform_;
  posture_transform_->SetMatrix(origin_matrix);

  double v[3] { p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2] };
  double normal[3];

  switch (direction) {
  case 0: {
    normal[0] = origin_matrix->GetElement(0, 0);
    normal[1] = origin_matrix->GetElement(1, 0);
    normal[2] = origin_matrix->GetElement(2, 0);
    std::cout << "normals: \n"; 
    for(int i = 0; i < 3; ++i){
      std::cout << normal[i] << " ";
    }
    std::cout << "\n";
    auto distance = vtkMath::Dot(v, normal);
    posture_transform_->Translate(distance, 0, 0);
  } break;
  case 1: {
    normal[0] = origin_matrix->GetElement(0, 1);
    normal[1] = origin_matrix->GetElement(1, 1);
    normal[2] = origin_matrix->GetElement(2, 1);
    std::cout << "normals: \n"; 
    for(int i = 0; i < 3; ++i){
      std::cout << normal[i] << " ";
    }
    std::cout << "\n";
    auto distance = vtkMath::Dot(v, normal);
    posture_transform_->Translate(0, distance, 0);
  } break;
  case 2: {
    normal[0] = origin_matrix->GetElement(0, 2);
    normal[1] = origin_matrix->GetElement(1, 2);
    normal[2] = origin_matrix->GetElement(2, 2);
    std::cout << "normals: \n"; 
    for(int i = 0; i < 3; ++i){
      std::cout << normal[i] << " ";
    }
    std::cout << "\n";
    auto distance = vtkMath::Dot(v, normal);
    posture_transform_->Translate(0, 0, distance);
  } break;
  default:
    return;
  }

  posture_transform_->Update();
  result_matrix->DeepCopy(posture_transform_->GetMatrix());
}

void vtkCustomRepresentation::GetPlaneProjection(double* a, double* b, double* c) {
  double sum = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
  c[0] = a[0] - sum*b[0];
  c[1] = a[1] - sum*b[1];
  c[2] = a[2] - sum*b[2];
  return;
}

void vtkCustomRepresentation::SetTargetActor(vtkActor** actor) {
  this->TargetActor = *actor;

  // matrix for rotate and move
  vtkNew<vtkMatrix4x4> Initial;
  this->TargetActor->SetUserMatrix(Initial);
}

void vtkCustomRepresentation::CreateDefaultProperties() {
  this->CenterProperty = vtkProperty::New();
  this->CenterProperty->SetColor(0.5, 0.5, 0.5);
  this->CenterProperty->SetAmbient(1.0);
  this->CenterProperty->SetOpacity(1.0);

  this->SelectedCenterProperty = vtkProperty::New();
  this->SelectedCenterProperty->SetAmbient(1.0);
  this->SelectedCenterProperty->SetDiffuse(1.0);
  this->SelectedCenterProperty->SetSpecular(0.0);
  this->SelectedCenterProperty->SetColor(1.0, 1.0, 0.0);
  this->SelectedCenterProperty->SetOpacity(0.5);
  // this->SelectedCenterProperty->SetOpacity(0.75);

  this->CircleProperty = new vtkProperty* [3];
  this->AxisProperty = new vtkProperty* [3];
  this->SelectedCircleProperty = new vtkProperty* [3];
  this->SelectedAxisProperty = new vtkProperty* [3];
  for(int i = 0; i < 3; ++i){
    this->CircleProperty[i] = vtkProperty::New();
    this->CircleProperty[i]->SetAmbient(1.0);
    // this->CircleProperty[i]->SetDiffuse(0.5);
    // this->CircleProperty[i]->SetSpecular(0.1);
    this->CircleProperty[i]->SetRepresentationToWireframe();
    this->CircleProperty[i]->SetLineWidth(this->AxisLineWidth);
    this->CircleProperty[i]->SetColor(this->AxisColor[i][0], this->AxisColor[i][1], this->AxisColor[i][2]);

    this->SelectedCircleProperty[i] = vtkProperty::New();
    this->SelectedCircleProperty[i]->SetAmbient(1.0);
    // this->SelectedCircleProperty[i]->SetDiffuse(1.0);
    // this->SelectedCircleProperty[i]->SetSpecular(1.0);
    this->SelectedCircleProperty[i]->SetColor(1, 1, 0);
    this->SelectedCircleProperty[i]->SetRepresentationToWireframe();
    this->SelectedCircleProperty[i]->SetLineWidth(this->AxisLineWidth);
    // this->SelectedCircleProperty[i]->SetColor(this->AxisColor[i][0], this->AxisColor[i][1], this->AxisColor[i][2]);

    this->AxisProperty[i] = vtkProperty::New();
    // this->LineProperty[i]->SetAmbient(1.0);
    this->AxisProperty[i]->SetLineWidth(this->AxisLineWidth);
    this->AxisProperty[i]->SetColor(this->AxisColor[i][0], this->AxisColor[i][1], this->AxisColor[i][2]);
  
    this->SelectedAxisProperty[i] = vtkProperty::New();
    this->SelectedAxisProperty[i]->SetAmbient(1.0);
    this->SelectedAxisProperty[i]->SetLineWidth(this->AxisLineWidth);
    // this->SelectedLineProperty[i]->EdgeVisibilityOn();
    // this->SelectedLineProperty[i]->SetColor(this->AxisColor[i][0], this->AxisColor[i][1], this->AxisColor[i][2]);
    this->SelectedAxisProperty[i]->SetColor(0.5, 0.5, 0.5);
  }
}

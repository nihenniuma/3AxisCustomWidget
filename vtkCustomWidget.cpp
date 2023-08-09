#include "vtkCustomWidget.h"
#include <vtkCommand.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkCallbackCommand.h>
#include <vtkWidgetCallbackMapper.h>
#include <vtkEvent.h>
#include <vtkWidgetEvent.h>
#include <vtkEventData.h>

vtkStandardNewMacro(vtkCustomWidget);

vtkCustomWidget::vtkCustomWidget() {
  // ----------------Need Test----------------------------
  this->CreateDefaultRepresentation(); 
  // ----------------Need Test----------------------------
  
  this->WidgetState = vtkCustomWidget::Start;

  this->ScalingEnabled = true;
  this->RotationEnabled = true;
  this->TranslationEnabled = true;
  this->CircleRotationEnabled  = true;
  this->LineTranslationEnabled  = true;

  this->CallbackMapper->SetCallbackMethod(vtkCommand::LeftButtonPressEvent,
                                        vtkEvent::NoModifier,
                                        0, 0, nullptr,
                                        vtkWidgetEvent::Select,
                                        this, vtkCustomWidget::SelectAction);

  this->CallbackMapper->SetCallbackMethod(vtkCommand::LeftButtonReleaseEvent,
                                        vtkEvent::NoModifier,
                                        0, 0, nullptr,
                                        vtkWidgetEvent::EndSelect,
                                        this, vtkCustomWidget::EndSelectAction);

  this->CallbackMapper->SetCallbackMethod(vtkCommand::MouseMoveEvent,
                                        vtkWidgetEvent::Rotate,
                                        this, vtkCustomWidget::RotateAction);

    this->CallbackMapper->SetCallbackMethod(vtkCommand::MouseMoveEvent,
                                        vtkWidgetEvent::Move,
                                        this, vtkCustomWidget::MoveAxisAction);

  // {
  //   vtkNew<vtkEventDataButton3D> ed;
  //   ed->SetDevice(vtkEventDataDevice::RightController);
  //   ed->SetInput(vtkEventDataDeviceInput::Trigger);
  //   ed->SetAction(vtkEventDataAction::Press);
  //   this->CallbackMapper->SetCallbackMethod(vtkCommand::Button3DEvent,
  //     ed.Get(), vtkWidgetEvent::Select3D,
  //     this, vtkCustomWidget::SelectAction3D);
  // }

  // {
  //   vtkNew<vtkEventDataButton3D> ed;
  //   ed->SetDevice(vtkEventDataDevice::RightController);
  //   ed->SetInput(vtkEventDataDeviceInput::Trigger);
  //   ed->SetAction(vtkEventDataAction::Release);
  //   this->CallbackMapper->SetCallbackMethod(vtkCommand::Button3DEvent,
  //     ed.Get(), vtkWidgetEvent::EndSelect3D,
  //     this, vtkCustomWidget::EndSelectAction3D);
  // }

  // {
  //   vtkNew<vtkEventDataMove3D> ed;
  //   ed->SetDevice(vtkEventDataDevice::RightController);
  //   this->CallbackMapper->SetCallbackMethod(vtkCommand::Move3DEvent,
  //     ed.Get(), vtkWidgetEvent::Rotate,
  //     this, vtkCustomWidget::RotateAction3D);
  // }
}

void vtkCustomWidget::SetTargetActor(vtkActor **actor) {
  if(this->WidgetRep){
    std::cout << "set target actor!\n"; 
    reinterpret_cast<vtkCustomRepresentation*>(this->WidgetRep)->SetTargetActor(actor);
  }
}

void vtkCustomWidget::CreateDefaultRepresentation() {
  if (!this->WidgetRep)
  {
    this->WidgetRep = vtkCustomRepresentation::New();
  }
}

void vtkCustomWidget::SelectAction(vtkAbstractWidget *w) {
  // We are in a static method, cast to ourself
  vtkCustomWidget *self = reinterpret_cast<vtkCustomWidget*>(w);

  // Get the event position
  int X = self->Interactor->GetEventPosition()[0];
  int Y = self->Interactor->GetEventPosition()[1];

  // Okay, make sure that the pick is in the current renderer
  if ( !self->CurrentRenderer ||
       !self->CurrentRenderer->IsInViewport(X,Y) )
  {
    self->WidgetState = vtkCustomWidget::Start;
    return;
  }

  // Begin the widget interaction which has the side effect of setting the
  // interaction state.
  double e[2];
  e[0] = static_cast<double>(X);
  e[1] = static_cast<double>(Y);
  self->WidgetRep->StartWidgetInteraction(e);
  int interactionState = self->WidgetRep->GetInteractionState();
  if ( interactionState == vtkCustomRepresentation::Outside )
  {
    return;
  }

  if((interactionState == vtkCustomRepresentation::RotateAroundX || 
      interactionState == vtkCustomRepresentation::RotateAroundY ||
      interactionState == vtkCustomRepresentation::RotateAroundZ)&&
      !self->RotationEnabled){
    return;
  }

  if((interactionState == vtkCustomRepresentation::MoveAlongX || 
      interactionState == vtkCustomRepresentation::MoveAlongZ ||
      interactionState == vtkCustomRepresentation::MoveAlongY ||
      interactionState == vtkCustomRepresentation::Moving)&&
      !self->TranslationEnabled){
    return;
  }

  self->WidgetState = vtkCustomWidget::Active;
  self->GrabFocus(self->EventCallbackCommand);

  // The SetInteractionState has the side effect of highlighting the widget
  reinterpret_cast<vtkCustomRepresentation*>(self->WidgetRep)->
    SetInteractionState(interactionState);

  // start the interaction
  self->EventCallbackCommand->SetAbortFlag(1);
  self->StartInteraction();
  self->InvokeEvent(vtkCommand::StartInteractionEvent,nullptr);
  self->Render();
}

void vtkCustomWidget::EndSelectAction(vtkAbstractWidget *w) {
  vtkCustomWidget *self = reinterpret_cast<vtkCustomWidget*>(w);
  if ( self->WidgetState == vtkCustomWidget::Start )
  {
    return;
  }

  // Return state to not active
  self->WidgetState = vtkCustomWidget::Start;
  reinterpret_cast<vtkCustomRepresentation*>(self->WidgetRep)->
    SetInteractionState(vtkCustomRepresentation::Outside);
  self->ReleaseFocus();

  self->EventCallbackCommand->SetAbortFlag(1);
  self->EndInteraction();
  self->InvokeEvent(vtkCommand::EndInteractionEvent,nullptr);
  self->Render();
}

void vtkCustomWidget::RotateAction(vtkAbstractWidget *w) {
  vtkCustomWidget *self = reinterpret_cast<vtkCustomWidget*>(w);

  // See whether we're active
  if ( self->WidgetState == vtkCustomWidget::Start )
  {
    return;
  }

  // compute some info we need for all cases
  int X = self->Interactor->GetEventPosition()[0];
  int Y = self->Interactor->GetEventPosition()[1];

  // Okay, adjust the representation
  double e[2];
  e[0] = static_cast<double>(X);
  e[1] = static_cast<double>(Y);
  self->WidgetRep->WidgetInteraction(e);

  // moving something
  self->EventCallbackCommand->SetAbortFlag(1);
  self->InvokeEvent(vtkCommand::InteractionEvent,nullptr);
  self->Render();
}

void vtkCustomWidget::MoveAxisAction(vtkAbstractWidget *w) {
  vtkCustomWidget *self = reinterpret_cast<vtkCustomWidget*>(w);

  // See whether we're active
  if ( self->WidgetState == vtkCustomWidget::Start )
  {
    return;
  }

  // compute some info we need for all cases
  int X = self->Interactor->GetEventPosition()[0];
  int Y = self->Interactor->GetEventPosition()[1];

  // Okay, adjust the representation
  double e[2];
  e[0] = static_cast<double>(X);
  e[1] = static_cast<double>(Y);
  self->WidgetRep->WidgetInteraction(e);

  // moving something
  self->EventCallbackCommand->SetAbortFlag(1);
  self->InvokeEvent(vtkCommand::InteractionEvent,nullptr);
  self->Render();
}

void vtkCustomWidget::RotateAction3D(vtkAbstractWidget *w) {
  vtkCustomWidget *self = reinterpret_cast<vtkCustomWidget*>(w);
  if ( self->WidgetState == vtkCustomWidget::Start )
  {
    return;
  }

  // Get the event position
  int X = self->Interactor->GetEventPosition()[0];
  int Y = self->Interactor->GetEventPosition()[1];

  // Okay, make sure that the pick is in the current renderer
  if ( !self->CurrentRenderer ||
       !self->CurrentRenderer->IsInViewport(X,Y) )
  {
    self->WidgetState = vtkCustomWidget::Start;
    return;
  }

  // Begin the widget interaction which has the side effect of setting the
  // interaction state.
  // double e[2];
  // e[0] = static_cast<double>(X);
  // e[1] = static_cast<double>(Y);
  // self->WidgetRep->StartWidgetInteraction(e);

  // self->WidgetRep->ComplexInteraction(
  //   self->Interactor, self, vtkWidgetEvent::Rotate, self->CallData);



  // int interactionState = self->WidgetRep->GetInteractionState();
  // switch (interactionState)
  // {
  //   case vtkCustomRepresentation::RotateAroundX:
  //     // reinterpret_cast<vtkCustomRepresentation*>(self->WidgetRep)->Rotate(X, Y, );
  //     self->WidgetState = vtkCustomWidget::Start;
  //     reinterpret_cast<vtkCustomRepresentation*>(self->WidgetRep)->ComplexInteraction(vtkCustomRepresentation::RotateAroundX);
  //     break;

  //   case vtkCustomRepresentation::RotateAroundY:
  //     // reinterpret_cast<vtkCustomRepresentation*>(self->WidgetRep)->Rotate(X, Y, );
  //     self->WidgetState = vtkCustomWidget::Start;
  //     reinterpret_cast<vtkCustomRepresentation*>(self->WidgetRep)->ComplexInteraction(vtkCustomRepresentation::RotateAroundY);
  //     break;
    
  //   case vtkCustomRepresentation::RotateAroundZ:
  //     // reinterpret_cast<vtkCustomRepresentation*>(self->WidgetRep)->Rotate(X, Y, );
  //     self->WidgetState = vtkCustomWidget::Start;
  //     reinterpret_cast<vtkCustomRepresentation*>(self->WidgetRep)->ComplexInteraction(vtkCustomRepresentation::RotateAroundZ);
  //     break;
    
  //   default:
  //     break;
  // }


  self->ReleaseFocus();

  self->EventCallbackCommand->SetAbortFlag(1);
  self->EndInteraction();
  self->InvokeEvent(vtkCommand::EndInteractionEvent,nullptr);
  self->Render();

}

void vtkCustomWidget::SelectAction3D(vtkAbstractWidget *w) {
  vtkCustomWidget *self = reinterpret_cast<vtkCustomWidget*>(w);

  // We want to compute an orthogonal vector to the plane that has been selected
  int interactionState = self->WidgetRep->ComputeComplexInteractionState(
    self->Interactor, self, vtkWidgetEvent::Select3D, self->CallData);

  if ( interactionState == vtkCustomRepresentation::Outside )
  {
    return;
  }

  // Test for states that involve face or handle picking here so
  // selection highlighting doesn't happen if that interaction is disabled.
  // Non-handle-grabbing transformations are tested in the "Action" methods.

  // // Rotation
  // if (interactionState == vtkCustomRepresentation::Rotating
  //      && self->RotationEnabled == 0)
  // {
  //   return;
  // }
  // // Face Movement
  // if ((interactionState == vtkCustomRepresentation::MoveF0 ||
  //      interactionState == vtkCustomRepresentation::MoveF1 ||
  //      interactionState == vtkCustomRepresentation::MoveF2 ||
  //      interactionState == vtkCustomRepresentation::MoveF3 ||
  //      interactionState == vtkCustomRepresentation::MoveF4 ||
  //      interactionState == vtkCustomRepresentation::MoveF5)
  //       && self->MoveFacesEnabled == 0)
  // {
  //   return;
  // }
  // // Translation
  // if (interactionState == vtkCustomRepresentation::Translating
  //      && self->TranslationEnabled == 0)
  // {
  //   return;
  // }

  // We are definitely selected
  if ( ! self->Parent )
  {
    self->GrabFocus(self->EventCallbackCommand);
  }

  self->WidgetState = vtkCustomWidget::Active;
  self->WidgetRep->StartComplexInteraction(
    self->Interactor, self, vtkWidgetEvent::Select3D, self->CallData);

  self->EventCallbackCommand->SetAbortFlag(1);
  self->StartInteraction();
  self->InvokeEvent(vtkCommand::StartInteractionEvent,nullptr);
}

void vtkCustomWidget::EndSelectAction3D(vtkAbstractWidget *w) {
  vtkCustomWidget *self = reinterpret_cast<vtkCustomWidget*>(w);

  if ( self->WidgetState != vtkCustomWidget::Active  ||
       self->WidgetRep->GetInteractionState() == vtkCustomRepresentation::Outside )
  {
    return;
  }

  // Return state to not selected
  self->WidgetRep->EndComplexInteraction(
    self->Interactor, self, vtkWidgetEvent::Select3D, self->CallData);

  self->WidgetState = vtkCustomWidget::Start;
  if ( ! self->Parent )
  {
    self->ReleaseFocus();
  }

  self->EventCallbackCommand->SetAbortFlag(1);
  self->EndInteraction();
  self->InvokeEvent(vtkCommand::EndInteractionEvent,nullptr);
}

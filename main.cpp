#include "vtkCustomWidget.h"
#include "vtkAxisWidget.h"

#include <vtkConeSource.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkBoxWidget2.h>
// #include <vtkBoxWidget.h>
// #include <vtkAffineWidget.h>

int main()
{
  vtkNew<vtkRenderer> renderer;
  vtkNew<vtkRenderWindow> renWin;
  vtkNew<vtkRenderWindowInteractor> renWinI;

  renderer->GetActiveCamera()->SetPosition(8, 8, 8);
  renderer->SetBackground(0.5, 0.5, 0.5);
  renWin->AddRenderer(renderer);
  renWinI->SetRenderWindow(renWin);
  renWinI->Initialize();
  // Widget的创建需要在 renWinI->SetRenderWindow(renWin) 和 renWinI->Initialize() 之后
  // 创建可移动坐标轴

  // vtkAxisWidget* axes = vtkAxisWidget::New();
  vtkCustomWidget* axes = vtkCustomWidget::New();
  // vtkBoxWidget2 *axes = vtkBoxWidget2::New();
  // vtkCustomRepresentation* re = vtkCustomRepresentation::New();
  // axes->SetRepresentation(re);
  axes->SetInteractor(renWinI);
  // axes->SetEnabled(1);

  vtkNew<vtkConeSource> ConeSource;
  vtkNew<vtkPolyDataMapper> ConeMapper;
  vtkActor* ConeActor = vtkActor::New();
  ConeMapper->SetInputConnection(ConeSource->GetOutputPort());
  ConeActor->SetMapper(ConeMapper);
  
  axes->SetTargetActor(&ConeActor);
  axes->On();

  // vtkBoxWidget2* box = vtkBoxWidget2::New();
  // box->SetInteractor(renWinI);
  // box->On();

  renWin->Render();
  renWinI->Start();

  return EXIT_SUCCESS;
}

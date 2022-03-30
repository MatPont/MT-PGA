#!/usr/bin/env python
import vtk
from vtk.util.misc import vtkGetDataRoot
VTK_DATA_ROOT = vtkGetDataRoot()

# Create the RenderWindow, Renderer and both Actors
#
ren1 = vtk.vtkRenderer()
renWin = vtk.vtkRenderWindow()
renWin.AddRenderer(ren1)
iren = vtk.vtkRenderWindowInteractor()
iren.SetRenderWindow(renWin)
# create a camera model
camCS = vtk.vtkConeSource()
camCS.SetHeight(1.5)
camCS.SetResolution(12)
camCS.SetRadius(0.4)
camCBS = vtk.vtkCubeSource()
camCBS.SetXLength(1.5)
camCBS.SetZLength(0.8)
camCBS.SetCenter(0.4,0,0)
camAPD = vtk.vtkAppendFilter()
camAPD.AddInputConnection(camCS.GetOutputPort())
camAPD.AddInputConnection(camCBS.GetOutputPort())
camMapper = vtk.vtkDataSetMapper()
camMapper.SetInputConnection(camAPD.GetOutputPort())
camActor = vtk.vtkLODActor()
camActor.SetMapper(camMapper)
camActor.SetScale(2,2,2)
# draw the arrows
pd = vtk.vtkPolyData()
ca = vtk.vtkCellArray()
fp = vtk.vtkPoints()
fp.InsertNextPoint(0,1,0)
fp.InsertNextPoint(8,1,0)
fp.InsertNextPoint(8,2,0)
fp.InsertNextPoint(10,0.01,0)
fp.InsertNextPoint(8,-2,0)
fp.InsertNextPoint(8,-1,0)
fp.InsertNextPoint(0,-1,0)
ca.InsertNextCell(7)
ca.InsertCellPoint(0)
ca.InsertCellPoint(1)
ca.InsertCellPoint(2)
ca.InsertCellPoint(3)
ca.InsertCellPoint(4)
ca.InsertCellPoint(5)
ca.InsertCellPoint(6)
pd.SetPoints(fp)
pd.SetPolys(ca)
pd2 = vtk.vtkPolyData()
ca2 = vtk.vtkCellArray()
fp2 = vtk.vtkPoints()
fp2.InsertNextPoint(0,1,0)
fp2.InsertNextPoint(8,1,0)
fp2.InsertNextPoint(8,2,0)
fp2.InsertNextPoint(10,0.01,0)
#prevents degenerate triangles
ca2.InsertNextCell(4)
ca2.InsertCellPoint(0)
ca2.InsertCellPoint(1)
ca2.InsertCellPoint(2)
ca2.InsertCellPoint(3)
pd2.SetPoints(fp2)
pd2.SetLines(ca2)
arrowIM = vtk.vtkImplicitModeller()
arrowIM.SetInputData(pd)
arrowIM.SetSampleDimensions(50,20,8)
arrowCF = vtk.vtkContourFilter()
arrowCF.SetInputConnection(arrowIM.GetOutputPort())
arrowCF.SetValue(0,0.2)
arrowWT = vtk.vtkWarpTo()
arrowWT.SetInputConnection(arrowCF.GetOutputPort())
arrowWT.SetPosition(5,0,5)
arrowWT.SetScaleFactor(0.85)
arrowWT.AbsoluteOn()
arrowT = vtk.vtkTransform()
arrowT.RotateY(60)
arrowT.Translate(-1.33198,0,-1.479)
arrowT.Scale(1,0.5,1)
arrowTF = vtk.vtkTransformFilter()
arrowTF.SetInputConnection(arrowWT.GetOutputPort())
arrowTF.SetTransform(arrowT)
arrowMapper = vtk.vtkDataSetMapper()
arrowMapper.SetInputConnection(arrowTF.GetOutputPort())
arrowMapper.ScalarVisibilityOff()
# draw the azimuth arrows
a1Actor = vtk.vtkLODActor()
a1Actor.SetMapper(arrowMapper)
a1Actor.RotateZ(180)
a1Actor.SetPosition(1,0,-1)
a1Actor.GetProperty().SetColor(1,0.3,0.3)
a1Actor.GetProperty().SetSpecularColor(1,1,1)
a1Actor.GetProperty().SetSpecular(0.3)
a1Actor.GetProperty().SetSpecularPower(20)
a1Actor.GetProperty().SetAmbient(0.2)
a1Actor.GetProperty().SetDiffuse(0.8)
a2Actor = vtk.vtkLODActor()
a2Actor.SetMapper(arrowMapper)
a2Actor.RotateZ(180)
a2Actor.RotateX(180)
a2Actor.SetPosition(1,0,1)
a2Actor.GetProperty().SetColor(1,0.3,0.3)
a2Actor.GetProperty().SetSpecularColor(1,1,1)
a2Actor.GetProperty().SetSpecular(0.3)
a2Actor.GetProperty().SetSpecularPower(20)
a2Actor.GetProperty().SetAmbient(0.2)
a2Actor.GetProperty().SetDiffuse(0.8)
# draw the elevation arrows
a3Actor = vtk.vtkLODActor()
a3Actor.SetMapper(arrowMapper)
a3Actor.RotateZ(180)
a3Actor.RotateX(90)
a3Actor.SetPosition(1,-1,0)
a3Actor.GetProperty().SetColor(0.3,1,0.3)
a3Actor.GetProperty().SetSpecularColor(1,1,1)
a3Actor.GetProperty().SetSpecular(0.3)
a3Actor.GetProperty().SetSpecularPower(20)
a3Actor.GetProperty().SetAmbient(0.2)
a3Actor.GetProperty().SetDiffuse(0.8)
a4Actor = vtk.vtkLODActor()
a4Actor.SetMapper(arrowMapper)
a4Actor.RotateZ(180)
a4Actor.RotateX(-90)
a4Actor.SetPosition(1,1,0)
a4Actor.GetProperty().SetColor(0.3,1,0.3)
a4Actor.GetProperty().SetSpecularColor(1,1,1)
a4Actor.GetProperty().SetSpecular(0.3)
a4Actor.GetProperty().SetSpecularPower(20)
a4Actor.GetProperty().SetAmbient(0.2)
a4Actor.GetProperty().SetDiffuse(0.8)
# draw the DOP
arrowT2 = vtk.vtkTransform()
arrowT2.Scale(1,0.6,1)
arrowT2.RotateY(90)
arrowTF2 = vtk.vtkTransformPolyDataFilter()
arrowTF2.SetInputData(pd2)
arrowTF2.SetTransform(arrowT2)
arrowREF = vtk.vtkRotationalExtrusionFilter()
arrowREF.SetInputConnection(arrowTF2.GetOutputPort())
arrowREF.CappingOff()
arrowREF.SetResolution(30)
spikeMapper = vtk.vtkPolyDataMapper()
spikeMapper.SetInputConnection(arrowREF.GetOutputPort())
a5Actor = vtk.vtkLODActor()
a5Actor.SetMapper(spikeMapper)
a5Actor.SetScale(.3,.3,.6)
a5Actor.RotateY(90)
a5Actor.SetPosition(-2,0,0)
a5Actor.GetProperty().SetColor(1,0.3,1)
a5Actor.GetProperty().SetAmbient(0.2)
a5Actor.GetProperty().SetDiffuse(0.8)
# focal point
fps = vtk.vtkSphereSource()
fps.SetRadius(0.5)
fpMapper = vtk.vtkPolyDataMapper()
fpMapper.SetInputConnection(fps.GetOutputPort())
fpActor = vtk.vtkLODActor()
fpActor.SetMapper(fpMapper)
fpActor.SetPosition(-9,0,0)
fpActor.GetProperty().SetSpecularColor(1,1,1)
fpActor.GetProperty().SetSpecular(0.3)
fpActor.GetProperty().SetAmbient(0.2)
fpActor.GetProperty().SetDiffuse(0.8)
fpActor.GetProperty().SetSpecularPower(20)
# create the roll arrows
arrowWT2 = vtk.vtkWarpTo()
arrowWT2.SetInputConnection(arrowCF.GetOutputPort())
arrowWT2.SetPosition(5,0,2.5)
arrowWT2.SetScaleFactor(0.95)
arrowWT2.AbsoluteOn()
arrowT3 = vtk.vtkTransform()
arrowT3.Translate(-2.50358,0,-1.70408)
arrowT3.Scale(0.5,0.3,1)
arrowTF3 = vtk.vtkTransformFilter()
arrowTF3.SetInputConnection(arrowWT2.GetOutputPort())
arrowTF3.SetTransform(arrowT3)
arrowMapper2 = vtk.vtkDataSetMapper()
arrowMapper2.SetInputConnection(arrowTF3.GetOutputPort())
arrowMapper2.ScalarVisibilityOff()
# draw the roll arrows
a6Actor = vtk.vtkLODActor()
a6Actor.SetMapper(arrowMapper2)
a6Actor.RotateZ(90)
a6Actor.SetPosition(-4,0,0)
a6Actor.SetScale(1.5,1.5,1.5)
a6Actor.GetProperty().SetColor(1,1,0.3)
a6Actor.GetProperty().SetSpecularColor(1,1,1)
a6Actor.GetProperty().SetSpecular(0.3)
a6Actor.GetProperty().SetSpecularPower(20)
a6Actor.GetProperty().SetAmbient(0.2)
a6Actor.GetProperty().SetDiffuse(0.8)
# Add the actors to the renderer, set the background and size
ren1.AddActor(camActor)
ren1.AddActor(a1Actor)
ren1.AddActor(a2Actor)
ren1.AddActor(a3Actor)
ren1.AddActor(a4Actor)
ren1.AddActor(a5Actor)
ren1.AddActor(a6Actor)
ren1.AddActor(fpActor)
ren1.SetBackground(0.1,0.2,0.4)
renWin.SetSize(300,300)
# render the image
ren1.ResetCamera()
cam1 = ren1.GetActiveCamera()
cam1.Zoom(1.5)
cam1.Azimuth(150)
cam1.Elevation(30)
iren.Initialize()
# prevent the tk window from showing up then start the event loop
# for testing
threshold = 15
# --- end of script --

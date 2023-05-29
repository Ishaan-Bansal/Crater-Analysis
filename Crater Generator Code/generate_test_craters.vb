' NX 1934
' Journal created by ibansal2 on Tue Mar 21 14:12:35 2023 Central Daylight Time
'
Imports System
Imports NXOpen

Module NXJournal
Sub Main (ByVal args() As String) 

For i As Integer = 1 to 6
    Dim radius As Double = i * 20.0
    Dim r As String = "" & radius
    For j As Integer = 1 to 7
        Dim depth As Double = j * 5.0
        Dim d As String = "" & depth

        Dim theSession As NXOpen.Session = NXOpen.Session.GetSession()
        Dim workPart As NXOpen.Part = theSession.Parts.Work

        Dim displayPart As NXOpen.Part = theSession.Parts.Display

        Dim markId1 As NXOpen.Session.UndoMarkId = Nothing
        markId1 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Visible, "Edit Sketch")

        theSession.BeginTaskEnvironment()

        Dim markId2 As NXOpen.Session.UndoMarkId = Nothing
        markId2 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Visible, "Enter Sketch")

        Dim sketchFeature1 As NXOpen.Features.SketchFeature = CType(workPart.Features.FindObject("SKETCH(1)"), NXOpen.Features.SketchFeature)

        Dim sketch1 As NXOpen.Sketch = CType(sketchFeature1.FindObject("SKETCH_000"), NXOpen.Sketch)

        sketch1.Activate(NXOpen.Sketch.ViewReorient.True)

        theSession.DeleteUndoMarksUpToMark(markId2, Nothing, True)

        Dim markId3 As NXOpen.Session.UndoMarkId = Nothing
        markId3 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Visible, "Open Sketch")

        Dim markId4 As NXOpen.Session.UndoMarkId = Nothing
        markId4 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Visible, "Edit Object Origin")

        Dim parallelDimension1 As NXOpen.Annotations.ParallelDimension = CType(theSession.ActiveSketch.FindObject("HANDLE R-2975"), NXOpen.Annotations.ParallelDimension)

        parallelDimension1.LeaderOrientation = NXOpen.Annotations.LeaderOrientation.FromRight

        parallelDimension1.IsOriginCentered = False

        Dim origin1 As NXOpen.Point3d = New NXOpen.Point3d(-16.380239962893782, 0.0, 9.6363950961766314)
        parallelDimension1.AnnotationOrigin = origin1

        Dim nErrs1 As Integer = Nothing
        nErrs1 = theSession.UpdateManager.DoUpdate(markId4)

        Dim markId5 As NXOpen.Session.UndoMarkId = Nothing
        markId5 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Visible, "Edit Dimension")

        Dim markId6 As NXOpen.Session.UndoMarkId = Nothing
        markId6 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Invisible, "Start")

        Dim sketchLinearDimensionBuilder1 As NXOpen.SketchLinearDimensionBuilder = Nothing
        sketchLinearDimensionBuilder1 = workPart.Sketches.CreateLinearDimensionBuilder(parallelDimension1)

        theSession.SetUndoMarkName(markId6, "Linear Dimension Dialog")

        theSession.SetUndoMarkVisibility(markId6, Nothing, NXOpen.Session.MarkVisibility.Visible)

        sketchLinearDimensionBuilder1.Origin.SetInferRelativeToGeometry(True)

        Dim dimensionlinearunits1 As NXOpen.Annotations.DimensionUnit = Nothing
        dimensionlinearunits1 = sketchLinearDimensionBuilder1.Style.UnitsStyle.DimensionLinearUnits

        Dim dimensionlinearunits2 As NXOpen.Annotations.DimensionUnit = Nothing
        dimensionlinearunits2 = sketchLinearDimensionBuilder1.Style.UnitsStyle.DimensionLinearUnits

        sketchLinearDimensionBuilder1.Style.OrdinateStyle.DoglegCreationOption = NXOpen.Annotations.OrdinateDoglegCreationOption.No

        Dim dimensionlinearunits3 As NXOpen.Annotations.DimensionUnit = Nothing
        dimensionlinearunits3 = sketchLinearDimensionBuilder1.Style.UnitsStyle.DimensionLinearUnits

        Dim dimensionlinearunits4 As NXOpen.Annotations.DimensionUnit = Nothing
        dimensionlinearunits4 = sketchLinearDimensionBuilder1.Style.UnitsStyle.DimensionLinearUnits

        Dim dimensionlinearunits5 As NXOpen.Annotations.DimensionUnit = Nothing
        dimensionlinearunits5 = sketchLinearDimensionBuilder1.Style.UnitsStyle.DimensionLinearUnits

        Dim dimensionlinearunits6 As NXOpen.Annotations.DimensionUnit = Nothing
        dimensionlinearunits6 = sketchLinearDimensionBuilder1.Style.UnitsStyle.DimensionLinearUnits

        sketchLinearDimensionBuilder1.Origin.SetInferRelativeToGeometry(True)

        sketchLinearDimensionBuilder1.Origin.SetInferRelativeToGeometry(True)

        Dim nullNXOpen_Direction As NXOpen.Direction = Nothing

        sketchLinearDimensionBuilder1.Measurement.Direction = nullNXOpen_Direction

        Dim nullNXOpen_View As NXOpen.View = Nothing

        sketchLinearDimensionBuilder1.Measurement.DirectionView = nullNXOpen_View

        sketchLinearDimensionBuilder1.Origin.SetInferRelativeToGeometry(True)

        Dim dimensionlinearunits7 As NXOpen.Annotations.DimensionUnit = Nothing
        dimensionlinearunits7 = sketchLinearDimensionBuilder1.Style.UnitsStyle.DimensionLinearUnits

        Dim dimensionlinearunits8 As NXOpen.Annotations.DimensionUnit = Nothing
        dimensionlinearunits8 = sketchLinearDimensionBuilder1.Style.UnitsStyle.DimensionLinearUnits

        Dim dimensionlinearunits9 As NXOpen.Annotations.DimensionUnit = Nothing
        dimensionlinearunits9 = sketchLinearDimensionBuilder1.Style.UnitsStyle.DimensionLinearUnits

        Dim dimensionlinearunits10 As NXOpen.Annotations.DimensionUnit = Nothing
        dimensionlinearunits10 = sketchLinearDimensionBuilder1.Style.UnitsStyle.DimensionLinearUnits

        Dim dimensionlinearunits11 As NXOpen.Annotations.DimensionUnit = Nothing
        dimensionlinearunits11 = sketchLinearDimensionBuilder1.Style.UnitsStyle.DimensionLinearUnits

        Dim dimensionlinearunits12 As NXOpen.Annotations.DimensionUnit = Nothing
        dimensionlinearunits12 = sketchLinearDimensionBuilder1.Style.UnitsStyle.DimensionLinearUnits

        Dim dimensionlinearunits13 As NXOpen.Annotations.DimensionUnit = Nothing
        dimensionlinearunits13 = sketchLinearDimensionBuilder1.Style.UnitsStyle.DimensionLinearUnits

        Dim dimensionlinearunits14 As NXOpen.Annotations.DimensionUnit = Nothing
        dimensionlinearunits14 = sketchLinearDimensionBuilder1.Style.UnitsStyle.DimensionLinearUnits

        Dim dimensionlinearunits15 As NXOpen.Annotations.DimensionUnit = Nothing
        dimensionlinearunits15 = sketchLinearDimensionBuilder1.Style.UnitsStyle.DimensionLinearUnits

        Dim dimensionlinearunits16 As NXOpen.Annotations.DimensionUnit = Nothing
        dimensionlinearunits16 = sketchLinearDimensionBuilder1.Style.UnitsStyle.DimensionLinearUnits

        Dim dimensionlinearunits17 As NXOpen.Annotations.DimensionUnit = Nothing
        dimensionlinearunits17 = sketchLinearDimensionBuilder1.Style.UnitsStyle.DimensionLinearUnits

        Dim dimensionlinearunits18 As NXOpen.Annotations.DimensionUnit = Nothing
        dimensionlinearunits18 = sketchLinearDimensionBuilder1.Style.UnitsStyle.DimensionLinearUnits

        Dim markId7 As NXOpen.Session.UndoMarkId = Nothing
        markId7 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Invisible, "Linear Dimension")

        theSession.DeleteUndoMark(markId7, Nothing)

        Dim markId8 As NXOpen.Session.UndoMarkId = Nothing
        markId8 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Invisible, "Linear Dimension")

        sketchLinearDimensionBuilder1.Driving.ExpressionValue.SetFormula(r)

        sketchLinearDimensionBuilder1.Driving.ExpressionMode = NXOpen.Annotations.DrivingValueBuilder.DrivingExpressionMode.KeepExpression

        sketchLinearDimensionBuilder1.Origin.SetInferRelativeToGeometry(True)

        Dim nXObject1 As NXOpen.NXObject = Nothing
        nXObject1 = sketchLinearDimensionBuilder1.Commit()

        Dim taggedObject1 As NXOpen.TaggedObject = Nothing
        taggedObject1 = sketchLinearDimensionBuilder1.SecondAssociativity.Value

        Dim line1 As NXOpen.Line = CType(taggedObject1, NXOpen.Line)

        Dim point1_1 As NXOpen.Point3d = New NXOpen.Point3d(0.0, 0.0, 0.0)
        Dim point2_1 As NXOpen.Point3d = New NXOpen.Point3d(0.0, 0.0, 0.0)
        sketchLinearDimensionBuilder1.FirstAssociativity.SetValue(NXOpen.InferSnapType.SnapType.Start, line1, nullNXOpen_View, point1_1, Nothing, nullNXOpen_View, point2_1)

        Dim point1_2 As NXOpen.Point3d = New NXOpen.Point3d(-radius, 0.0, 0.0)
        Dim point2_2 As NXOpen.Point3d = New NXOpen.Point3d(0.0, 0.0, 0.0)
        sketchLinearDimensionBuilder1.SecondAssociativity.SetValue(NXOpen.InferSnapType.SnapType.End, line1, nullNXOpen_View, point1_2, Nothing, nullNXOpen_View, point2_2)

        sketchLinearDimensionBuilder1.Driving.ExpressionValue.SetFormula(r)

        theSession.SetUndoMarkName(markId8, "Linear Dimension - =")

        theSession.SetUndoMarkVisibility(markId8, Nothing, NXOpen.Session.MarkVisibility.Visible)

        theSession.SetUndoMarkVisibility(markId6, Nothing, NXOpen.Session.MarkVisibility.Invisible)

        Dim markId9 As NXOpen.Session.UndoMarkId = Nothing
        markId9 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Invisible, "Linear Dimension")

        Dim markId10 As NXOpen.Session.UndoMarkId = Nothing
        markId10 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Invisible, "Linear Dimension")

        Dim nXObject2 As NXOpen.NXObject = Nothing
        nXObject2 = sketchLinearDimensionBuilder1.Commit()

        theSession.DeleteUndoMark(markId10, Nothing)

        theSession.SetUndoMarkName(markId6, "Linear Dimension")

        Dim expression1 As NXOpen.Expression = sketchLinearDimensionBuilder1.Driving.ExpressionValue

        sketchLinearDimensionBuilder1.Destroy()

        theSession.DeleteUndoMark(markId9, Nothing)

        theSession.SetUndoMarkVisibility(markId6, Nothing, NXOpen.Session.MarkVisibility.Visible)

        theSession.DeleteUndoMark(markId8, Nothing)

        theSession.DeleteUndoMark(markId6, Nothing)

        Dim markId11 As NXOpen.Session.UndoMarkId = Nothing
        markId11 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Visible, "Edit Dimension")

        Dim markId12 As NXOpen.Session.UndoMarkId = Nothing
        markId12 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Invisible, "Start")

        Dim parallelDimension2 As NXOpen.Annotations.ParallelDimension = CType(theSession.ActiveSketch.FindObject("HANDLE R-3264"), NXOpen.Annotations.ParallelDimension)

        Dim sketchLinearDimensionBuilder2 As NXOpen.SketchLinearDimensionBuilder = Nothing
        sketchLinearDimensionBuilder2 = workPart.Sketches.CreateLinearDimensionBuilder(parallelDimension2)

        theSession.SetUndoMarkName(markId12, "Linear Dimension Dialog")

        theSession.SetUndoMarkVisibility(markId12, Nothing, NXOpen.Session.MarkVisibility.Visible)

        sketchLinearDimensionBuilder2.Origin.SetInferRelativeToGeometry(True)

        Dim dimensionlinearunits19 As NXOpen.Annotations.DimensionUnit = Nothing
        dimensionlinearunits19 = sketchLinearDimensionBuilder2.Style.UnitsStyle.DimensionLinearUnits

        Dim dimensionlinearunits20 As NXOpen.Annotations.DimensionUnit = Nothing
        dimensionlinearunits20 = sketchLinearDimensionBuilder2.Style.UnitsStyle.DimensionLinearUnits

        sketchLinearDimensionBuilder2.Style.OrdinateStyle.DoglegCreationOption = NXOpen.Annotations.OrdinateDoglegCreationOption.No

        Dim dimensionlinearunits21 As NXOpen.Annotations.DimensionUnit = Nothing
        dimensionlinearunits21 = sketchLinearDimensionBuilder2.Style.UnitsStyle.DimensionLinearUnits

        Dim dimensionlinearunits22 As NXOpen.Annotations.DimensionUnit = Nothing
        dimensionlinearunits22 = sketchLinearDimensionBuilder2.Style.UnitsStyle.DimensionLinearUnits

        Dim dimensionlinearunits23 As NXOpen.Annotations.DimensionUnit = Nothing
        dimensionlinearunits23 = sketchLinearDimensionBuilder2.Style.UnitsStyle.DimensionLinearUnits

        Dim dimensionlinearunits24 As NXOpen.Annotations.DimensionUnit = Nothing
        dimensionlinearunits24 = sketchLinearDimensionBuilder2.Style.UnitsStyle.DimensionLinearUnits

        sketchLinearDimensionBuilder2.Origin.SetInferRelativeToGeometry(True)

        sketchLinearDimensionBuilder2.Origin.SetInferRelativeToGeometry(True)

        sketchLinearDimensionBuilder2.Measurement.Direction = nullNXOpen_Direction

        sketchLinearDimensionBuilder2.Measurement.DirectionView = nullNXOpen_View

        sketchLinearDimensionBuilder2.Origin.SetInferRelativeToGeometry(True)

        Dim dimensionlinearunits25 As NXOpen.Annotations.DimensionUnit = Nothing
        dimensionlinearunits25 = sketchLinearDimensionBuilder2.Style.UnitsStyle.DimensionLinearUnits

        Dim dimensionlinearunits26 As NXOpen.Annotations.DimensionUnit = Nothing
        dimensionlinearunits26 = sketchLinearDimensionBuilder2.Style.UnitsStyle.DimensionLinearUnits

        Dim dimensionlinearunits27 As NXOpen.Annotations.DimensionUnit = Nothing
        dimensionlinearunits27 = sketchLinearDimensionBuilder2.Style.UnitsStyle.DimensionLinearUnits

        Dim dimensionlinearunits28 As NXOpen.Annotations.DimensionUnit = Nothing
        dimensionlinearunits28 = sketchLinearDimensionBuilder2.Style.UnitsStyle.DimensionLinearUnits

        Dim dimensionlinearunits29 As NXOpen.Annotations.DimensionUnit = Nothing
        dimensionlinearunits29 = sketchLinearDimensionBuilder2.Style.UnitsStyle.DimensionLinearUnits

        Dim dimensionlinearunits30 As NXOpen.Annotations.DimensionUnit = Nothing
        dimensionlinearunits30 = sketchLinearDimensionBuilder2.Style.UnitsStyle.DimensionLinearUnits

        Dim dimensionlinearunits31 As NXOpen.Annotations.DimensionUnit = Nothing
        dimensionlinearunits31 = sketchLinearDimensionBuilder2.Style.UnitsStyle.DimensionLinearUnits

        Dim dimensionlinearunits32 As NXOpen.Annotations.DimensionUnit = Nothing
        dimensionlinearunits32 = sketchLinearDimensionBuilder2.Style.UnitsStyle.DimensionLinearUnits

        Dim dimensionlinearunits33 As NXOpen.Annotations.DimensionUnit = Nothing
        dimensionlinearunits33 = sketchLinearDimensionBuilder2.Style.UnitsStyle.DimensionLinearUnits

        Dim dimensionlinearunits34 As NXOpen.Annotations.DimensionUnit = Nothing
        dimensionlinearunits34 = sketchLinearDimensionBuilder2.Style.UnitsStyle.DimensionLinearUnits

        Dim dimensionlinearunits35 As NXOpen.Annotations.DimensionUnit = Nothing
        dimensionlinearunits35 = sketchLinearDimensionBuilder2.Style.UnitsStyle.DimensionLinearUnits

        Dim dimensionlinearunits36 As NXOpen.Annotations.DimensionUnit = Nothing
        dimensionlinearunits36 = sketchLinearDimensionBuilder2.Style.UnitsStyle.DimensionLinearUnits

        ' ----------------------------------------------
        '   Dialog Begin Linear Dimension
        ' ----------------------------------------------
        Dim markId13 As NXOpen.Session.UndoMarkId = Nothing
        markId13 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Invisible, "Linear Dimension")

        theSession.DeleteUndoMark(markId13, Nothing)

        Dim markId14 As NXOpen.Session.UndoMarkId = Nothing
        markId14 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Invisible, "Linear Dimension")

        sketchLinearDimensionBuilder2.Driving.ExpressionValue.SetFormula(depth)

        sketchLinearDimensionBuilder2.Driving.ExpressionMode = NXOpen.Annotations.DrivingValueBuilder.DrivingExpressionMode.KeepExpression

        sketchLinearDimensionBuilder2.Origin.SetInferRelativeToGeometry(True)

        Dim nXObject3 As NXOpen.NXObject = Nothing
        nXObject3 = sketchLinearDimensionBuilder2.Commit()

        Dim taggedObject2 As NXOpen.TaggedObject = Nothing
        taggedObject2 = sketchLinearDimensionBuilder2.SecondAssociativity.Value

        Dim line2 As NXOpen.Line = CType(taggedObject2, NXOpen.Line)

        Dim point1_3 As NXOpen.Point3d = New NXOpen.Point3d(0.0, 0.0, 0.0)
        Dim point2_3 As NXOpen.Point3d = New NXOpen.Point3d(0.0, 0.0, 0.0)
        sketchLinearDimensionBuilder2.FirstAssociativity.SetValue(NXOpen.InferSnapType.SnapType.Start, line2, nullNXOpen_View, point1_3, Nothing, nullNXOpen_View, point2_3)

        Dim point1_4 As NXOpen.Point3d = New NXOpen.Point3d(0.0, 0.0, -depth)
        Dim point2_4 As NXOpen.Point3d = New NXOpen.Point3d(0.0, 0.0, 0.0)
        sketchLinearDimensionBuilder2.SecondAssociativity.SetValue(NXOpen.InferSnapType.SnapType.End, line2, nullNXOpen_View, point1_4, Nothing, nullNXOpen_View, point2_4)

        sketchLinearDimensionBuilder2.Driving.ExpressionValue.SetFormula(d)

        theSession.SetUndoMarkName(markId14, "Linear Dimension - =")

        theSession.SetUndoMarkVisibility(markId14, Nothing, NXOpen.Session.MarkVisibility.Visible)

        theSession.SetUndoMarkVisibility(markId12, Nothing, NXOpen.Session.MarkVisibility.Invisible)

        Dim markId15 As NXOpen.Session.UndoMarkId = Nothing
        markId15 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Invisible, "Linear Dimension")

        Dim markId16 As NXOpen.Session.UndoMarkId = Nothing
        markId16 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Invisible, "Linear Dimension")

        Dim nXObject4 As NXOpen.NXObject = Nothing
        nXObject4 = sketchLinearDimensionBuilder2.Commit()

        theSession.DeleteUndoMark(markId16, Nothing)

        theSession.SetUndoMarkName(markId12, "Linear Dimension")

        Dim expression2 As NXOpen.Expression = sketchLinearDimensionBuilder2.Driving.ExpressionValue

        sketchLinearDimensionBuilder2.Destroy()

        theSession.DeleteUndoMark(markId15, Nothing)

        theSession.SetUndoMarkVisibility(markId12, Nothing, NXOpen.Session.MarkVisibility.Visible)

        theSession.DeleteUndoMark(markId14, Nothing)

        theSession.DeleteUndoMark(markId12, Nothing)

        ' ----------------------------------------------
        '   Menu: Task->Finish Sketch
        ' ----------------------------------------------
        Dim revolve1 As NXOpen.Features.Revolve = CType(workPart.Features.FindObject("REVOLVED(2)"), NXOpen.Features.Revolve)

        Dim section1() As NXOpen.Section
        section1 = revolve1.GetSections()

        section1(0).CleanMappingData()

        theSession.Preferences.Sketch.SectionView = False

        Dim markId17 As NXOpen.Session.UndoMarkId = Nothing
        markId17 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Invisible, "Deactivate Sketch")

        theSession.ActiveSketch.Deactivate(NXOpen.Sketch.ViewReorient.True, NXOpen.Sketch.UpdateLevel.Model)

        theSession.DeleteUndoMarksSetInTaskEnvironment()

        theSession.EndTaskEnvironment()

        ' ----------------------------------------------
        '   Menu: File->Export->STL...
        ' ----------------------------------------------
        Dim markId18 As NXOpen.Session.UndoMarkId = Nothing
        markId18 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Visible, "Start")

        Dim sTLCreator1 As NXOpen.STLCreator = Nothing
        sTLCreator1 = theSession.DexManager.CreateStlCreator()

        sTLCreator1.AutoNormalGen = True

        sTLCreator1.ChordalTol = 0.080000000000000002

        sTLCreator1.AdjacencyTol = 0.080000000000000002

        theSession.SetUndoMarkName(markId18, "STL Export Dialog")

        Dim body1 As NXOpen.Body = CType(workPart.Bodies.FindObject("REVOLVED(2)"), NXOpen.Body)

        Dim added1 As Boolean = Nothing
        added1 = sTLCreator1.ExportSelectionBlock.Add(body1)

        Dim markId19 As NXOpen.Session.UndoMarkId = Nothing
        markId19 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Invisible, "STL Export")

        theSession.DeleteUndoMark(markId19, Nothing)

        Dim markId20 As NXOpen.Session.UndoMarkId = Nothing
        markId20 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Invisible, "STL Export")

        sTLCreator1.OutputFile = "\\ad.uillinois.edu\engr-ews\ibansal2\Desktop\Plumes_Research\Test Craters STLs\r" + r + "_d" + d + ".stl"

        Dim nXObject5 As NXOpen.NXObject = Nothing
        nXObject5 = sTLCreator1.Commit()

        theSession.DeleteUndoMark(markId20, Nothing)

        theSession.SetUndoMarkName(markId18, "STL Export")

        sTLCreator1.Destroy()
    Next
Next

' ----------------------------------------------
'   Menu: Tools->Journal->Stop Recording
' ----------------------------------------------

End Sub
End Module
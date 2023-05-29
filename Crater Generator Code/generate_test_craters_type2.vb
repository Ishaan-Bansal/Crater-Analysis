' NX 1934
' Journal created by ibansal2 on Fri Apr  7 15:05:34 2023 Central Daylight Time
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

        For k as Integer = 1 to 5
            Dim curvature as Double = k * 5.0
            Dim c As String = "" & curvature

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

            Dim radiusDimension1 As NXOpen.Annotations.RadiusDimension = CType(theSession.ActiveSketch.FindObject("Dimension p5"), NXOpen.Annotations.RadiusDimension)

            radiusDimension1.LeaderOrientation = NXOpen.Annotations.LeaderOrientation.FromLeft

            radiusDimension1.IsOriginCentered = False

            Dim origin1 As NXOpen.Point3d = New NXOpen.Point3d(-126.59722624013595, 0.0, -6.7617932190186272)
            radiusDimension1.AnnotationOrigin = origin1

            Dim nErrs1 As Integer = Nothing
            nErrs1 = theSession.UpdateManager.DoUpdate(markId4)

            Dim scaleAboutPoint1 As NXOpen.Point3d = New NXOpen.Point3d(11.444416866172528, -6.3700056141903536, 0.0)
            Dim viewCenter1 As NXOpen.Point3d = New NXOpen.Point3d(-11.444416866172528, 6.3700056141903536, 0.0)
            workPart.ModelingViews.WorkView.ZoomAboutPoint(1.25, scaleAboutPoint1, viewCenter1)

            Dim markId5 As NXOpen.Session.UndoMarkId = Nothing
            markId5 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Visible, "Edit Dimension")

            Dim markId6 As NXOpen.Session.UndoMarkId = Nothing
            markId6 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Invisible, "Start")

            Dim sketchRadialDimensionBuilder1 As NXOpen.SketchRadialDimensionBuilder = Nothing
            sketchRadialDimensionBuilder1 = workPart.Sketches.CreateRadialDimensionBuilder(radiusDimension1)

            theSession.SetUndoMarkName(markId6, "Radial Dimension Dialog")

            theSession.SetUndoMarkVisibility(markId6, Nothing, NXOpen.Session.MarkVisibility.Visible)

            sketchRadialDimensionBuilder1.Origin.SetInferRelativeToGeometry(True)

            Dim dimensionlinearunits1 As NXOpen.Annotations.DimensionUnit = Nothing
            dimensionlinearunits1 = sketchRadialDimensionBuilder1.Style.UnitsStyle.DimensionLinearUnits

            Dim dimensionlinearunits2 As NXOpen.Annotations.DimensionUnit = Nothing
            dimensionlinearunits2 = sketchRadialDimensionBuilder1.Style.UnitsStyle.DimensionLinearUnits

            sketchRadialDimensionBuilder1.Style.OrdinateStyle.DoglegCreationOption = NXOpen.Annotations.OrdinateDoglegCreationOption.No

            Dim dimensionlinearunits3 As NXOpen.Annotations.DimensionUnit = Nothing
            dimensionlinearunits3 = sketchRadialDimensionBuilder1.Style.UnitsStyle.DimensionLinearUnits

            Dim dimensionlinearunits4 As NXOpen.Annotations.DimensionUnit = Nothing
            dimensionlinearunits4 = sketchRadialDimensionBuilder1.Style.UnitsStyle.DimensionLinearUnits

            Dim dimensionlinearunits5 As NXOpen.Annotations.DimensionUnit = Nothing
            dimensionlinearunits5 = sketchRadialDimensionBuilder1.Style.UnitsStyle.DimensionLinearUnits

            Dim dimensionlinearunits6 As NXOpen.Annotations.DimensionUnit = Nothing
            dimensionlinearunits6 = sketchRadialDimensionBuilder1.Style.UnitsStyle.DimensionLinearUnits

            sketchRadialDimensionBuilder1.Origin.SetInferRelativeToGeometry(True)

            sketchRadialDimensionBuilder1.Origin.SetInferRelativeToGeometry(True)

            Dim nullNXOpen_Direction As NXOpen.Direction = Nothing

            sketchRadialDimensionBuilder1.Measurement.Direction = nullNXOpen_Direction

            Dim nullNXOpen_View As NXOpen.View = Nothing

            sketchRadialDimensionBuilder1.Measurement.DirectionView = nullNXOpen_View

            Dim dimensionlinearunits7 As NXOpen.Annotations.DimensionUnit = Nothing
            dimensionlinearunits7 = sketchRadialDimensionBuilder1.Style.UnitsStyle.DimensionLinearUnits

            Dim dimensionlinearunits8 As NXOpen.Annotations.DimensionUnit = Nothing
            dimensionlinearunits8 = sketchRadialDimensionBuilder1.Style.UnitsStyle.DimensionLinearUnits

            Dim dimensionlinearunits9 As NXOpen.Annotations.DimensionUnit = Nothing
            dimensionlinearunits9 = sketchRadialDimensionBuilder1.Style.UnitsStyle.DimensionLinearUnits

            Dim dimensionlinearunits10 As NXOpen.Annotations.DimensionUnit = Nothing
            dimensionlinearunits10 = sketchRadialDimensionBuilder1.Style.UnitsStyle.DimensionLinearUnits

            Dim dimensionlinearunits11 As NXOpen.Annotations.DimensionUnit = Nothing
            dimensionlinearunits11 = sketchRadialDimensionBuilder1.Style.UnitsStyle.DimensionLinearUnits

            Dim dimensionlinearunits12 As NXOpen.Annotations.DimensionUnit = Nothing
            dimensionlinearunits12 = sketchRadialDimensionBuilder1.Style.UnitsStyle.DimensionLinearUnits

            Dim dimensionlinearunits13 As NXOpen.Annotations.DimensionUnit = Nothing
            dimensionlinearunits13 = sketchRadialDimensionBuilder1.Style.UnitsStyle.DimensionLinearUnits

            Dim dimensionlinearunits14 As NXOpen.Annotations.DimensionUnit = Nothing
            dimensionlinearunits14 = sketchRadialDimensionBuilder1.Style.UnitsStyle.DimensionLinearUnits

            Dim dimensionlinearunits15 As NXOpen.Annotations.DimensionUnit = Nothing
            dimensionlinearunits15 = sketchRadialDimensionBuilder1.Style.UnitsStyle.DimensionLinearUnits

            Dim dimensionlinearunits16 As NXOpen.Annotations.DimensionUnit = Nothing
            dimensionlinearunits16 = sketchRadialDimensionBuilder1.Style.UnitsStyle.DimensionLinearUnits

            Dim dimensionlinearunits17 As NXOpen.Annotations.DimensionUnit = Nothing
            dimensionlinearunits17 = sketchRadialDimensionBuilder1.Style.UnitsStyle.DimensionLinearUnits

            Dim dimensionlinearunits18 As NXOpen.Annotations.DimensionUnit = Nothing
            dimensionlinearunits18 = sketchRadialDimensionBuilder1.Style.UnitsStyle.DimensionLinearUnits

            ' ----------------------------------------------
            '   Dialog Begin Radial Dimension
            ' ----------------------------------------------
            Dim markId7 As NXOpen.Session.UndoMarkId = Nothing
            markId7 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Invisible, "Radial Dimension")

            theSession.DeleteUndoMark(markId7, Nothing)

            Dim markId8 As NXOpen.Session.UndoMarkId = Nothing
            markId8 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Invisible, "Radial Dimension")

            sketchRadialDimensionBuilder1.Driving.ExpressionValue.SetFormula(c)

            sketchRadialDimensionBuilder1.Driving.ExpressionMode = NXOpen.Annotations.DrivingValueBuilder.DrivingExpressionMode.KeepExpression

            Dim nXObject1 As NXOpen.NXObject = Nothing
            nXObject1 = sketchRadialDimensionBuilder1.Commit()

            Dim taggedObject1 As NXOpen.TaggedObject = Nothing
            taggedObject1 = sketchRadialDimensionBuilder1.FirstAssociativity.Value

            Dim arc1 As NXOpen.Arc = CType(taggedObject1, NXOpen.Arc)

            Dim point1_1 As NXOpen.Point3d = New NXOpen.Point3d(-127.44295466168791, 0.0, -c)
            Dim point2_1 As NXOpen.Point3d = New NXOpen.Point3d(0.0, 0.0, 0.0)
            sketchRadialDimensionBuilder1.FirstAssociativity.SetValue(NXOpen.InferSnapType.SnapType.None, arc1, nullNXOpen_View, point1_1, Nothing, nullNXOpen_View, point2_1)

            sketchRadialDimensionBuilder1.Driving.ExpressionValue.SetFormula("10")

            theSession.SetUndoMarkName(markId8, "Radial Dimension - =")

            theSession.SetUndoMarkVisibility(markId8, Nothing, NXOpen.Session.MarkVisibility.Visible)

            theSession.SetUndoMarkVisibility(markId6, Nothing, NXOpen.Session.MarkVisibility.Invisible)

            Dim markId9 As NXOpen.Session.UndoMarkId = Nothing
            markId9 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Invisible, "Radial Dimension")

            Dim markId10 As NXOpen.Session.UndoMarkId = Nothing
            markId10 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Invisible, "Radial Dimension")

            Dim nXObject2 As NXOpen.NXObject = Nothing
            nXObject2 = sketchRadialDimensionBuilder1.Commit()

            theSession.DeleteUndoMark(markId10, Nothing)

            theSession.SetUndoMarkName(markId6, "Radial Dimension")

            Dim expression1 As NXOpen.Expression = sketchRadialDimensionBuilder1.Driving.ExpressionValue

            sketchRadialDimensionBuilder1.Destroy()

            theSession.DeleteUndoMark(markId9, Nothing)

            theSession.SetUndoMarkVisibility(markId6, Nothing, NXOpen.Session.MarkVisibility.Visible)

            theSession.DeleteUndoMark(markId8, Nothing)

            theSession.DeleteUndoMark(markId6, Nothing)

            Dim scaleAboutPoint2 As NXOpen.Point3d = New NXOpen.Point3d(4.3186478740274135, 3.8867830866246309, 0.0)
            Dim viewCenter2 As NXOpen.Point3d = New NXOpen.Point3d(-4.3186478740273841, -3.8867830866246309, 0.0)
            workPart.ModelingViews.WorkView.ZoomAboutPoint(0.80000000000000004, scaleAboutPoint2, viewCenter2)

            Dim scaleAboutPoint3 As NXOpen.Point3d = New NXOpen.Point3d(5.3983098425342311, 4.210681677176682, 0.0)
            Dim viewCenter3 As NXOpen.Point3d = New NXOpen.Point3d(-5.3983098425342311, -4.210681677176682, 0.0)
            workPart.ModelingViews.WorkView.ZoomAboutPoint(0.80000000000000004, scaleAboutPoint3, viewCenter3)

            Dim scaleAboutPoint4 As NXOpen.Point3d = New NXOpen.Point3d(6.7478873031678122, 4.453605620090741, 0.0)
            Dim viewCenter4 As NXOpen.Point3d = New NXOpen.Point3d(-6.7478873031678122, -4.453605620090741, 0.0)
            workPart.ModelingViews.WorkView.ZoomAboutPoint(0.80000000000000004, scaleAboutPoint4, viewCenter4)

            Dim scaleAboutPoint5 As NXOpen.Point3d = New NXOpen.Point3d(91.096478592764839, 40.993415366744138, 0.0)
            Dim viewCenter5 As NXOpen.Point3d = New NXOpen.Point3d(-91.096478592764839, -40.993415366744195, 0.0)
            workPart.ModelingViews.WorkView.ZoomAboutPoint(1.25, scaleAboutPoint5, viewCenter5)

            Dim scaleAboutPoint6 As NXOpen.Point3d = New NXOpen.Point3d(72.87718287421194, 32.794732293395313, 0.0)
            Dim viewCenter6 As NXOpen.Point3d = New NXOpen.Point3d(-72.877182874211883, -32.794732293395363, 0.0)
            workPart.ModelingViews.WorkView.ZoomAboutPoint(1.25, scaleAboutPoint6, viewCenter6)

            Dim scaleAboutPoint7 As NXOpen.Point3d = New NXOpen.Point3d(58.301746299369547, 26.235785834716214, 0.0)
            Dim viewCenter7 As NXOpen.Point3d = New NXOpen.Point3d(-58.301746299369434, -26.235785834716307, 0.0)
            workPart.ModelingViews.WorkView.ZoomAboutPoint(1.25, scaleAboutPoint7, viewCenter7)

            Dim markId11 As NXOpen.Session.UndoMarkId = Nothing
            markId11 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Visible, "Edit Dimension")

            Dim markId12 As NXOpen.Session.UndoMarkId = Nothing
            markId12 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Invisible, "Start")

            Dim parallelDimension1 As NXOpen.Annotations.ParallelDimension = CType(theSession.ActiveSketch.FindObject("HANDLE R-2975"), NXOpen.Annotations.ParallelDimension)

            Dim sketchLinearDimensionBuilder1 As NXOpen.SketchLinearDimensionBuilder = Nothing
            sketchLinearDimensionBuilder1 = workPart.Sketches.CreateLinearDimensionBuilder(parallelDimension1)

            theSession.SetUndoMarkName(markId12, "Linear Dimension Dialog")

            theSession.SetUndoMarkVisibility(markId12, Nothing, NXOpen.Session.MarkVisibility.Visible)

            sketchLinearDimensionBuilder1.Origin.SetInferRelativeToGeometry(True)

            Dim dimensionlinearunits19 As NXOpen.Annotations.DimensionUnit = Nothing
            dimensionlinearunits19 = sketchLinearDimensionBuilder1.Style.UnitsStyle.DimensionLinearUnits

            Dim dimensionlinearunits20 As NXOpen.Annotations.DimensionUnit = Nothing
            dimensionlinearunits20 = sketchLinearDimensionBuilder1.Style.UnitsStyle.DimensionLinearUnits

            sketchLinearDimensionBuilder1.Style.OrdinateStyle.DoglegCreationOption = NXOpen.Annotations.OrdinateDoglegCreationOption.No

            Dim dimensionlinearunits21 As NXOpen.Annotations.DimensionUnit = Nothing
            dimensionlinearunits21 = sketchLinearDimensionBuilder1.Style.UnitsStyle.DimensionLinearUnits

            Dim dimensionlinearunits22 As NXOpen.Annotations.DimensionUnit = Nothing
            dimensionlinearunits22 = sketchLinearDimensionBuilder1.Style.UnitsStyle.DimensionLinearUnits

            Dim dimensionlinearunits23 As NXOpen.Annotations.DimensionUnit = Nothing
            dimensionlinearunits23 = sketchLinearDimensionBuilder1.Style.UnitsStyle.DimensionLinearUnits

            Dim dimensionlinearunits24 As NXOpen.Annotations.DimensionUnit = Nothing
            dimensionlinearunits24 = sketchLinearDimensionBuilder1.Style.UnitsStyle.DimensionLinearUnits

            sketchLinearDimensionBuilder1.Origin.SetInferRelativeToGeometry(True)

            sketchLinearDimensionBuilder1.Origin.SetInferRelativeToGeometry(True)

            sketchLinearDimensionBuilder1.Measurement.Direction = nullNXOpen_Direction

            sketchLinearDimensionBuilder1.Measurement.DirectionView = nullNXOpen_View

            sketchLinearDimensionBuilder1.Origin.SetInferRelativeToGeometry(True)

            Dim dimensionlinearunits25 As NXOpen.Annotations.DimensionUnit = Nothing
            dimensionlinearunits25 = sketchLinearDimensionBuilder1.Style.UnitsStyle.DimensionLinearUnits

            Dim dimensionlinearunits26 As NXOpen.Annotations.DimensionUnit = Nothing
            dimensionlinearunits26 = sketchLinearDimensionBuilder1.Style.UnitsStyle.DimensionLinearUnits

            Dim dimensionlinearunits27 As NXOpen.Annotations.DimensionUnit = Nothing
            dimensionlinearunits27 = sketchLinearDimensionBuilder1.Style.UnitsStyle.DimensionLinearUnits

            Dim dimensionlinearunits28 As NXOpen.Annotations.DimensionUnit = Nothing
            dimensionlinearunits28 = sketchLinearDimensionBuilder1.Style.UnitsStyle.DimensionLinearUnits

            Dim dimensionlinearunits29 As NXOpen.Annotations.DimensionUnit = Nothing
            dimensionlinearunits29 = sketchLinearDimensionBuilder1.Style.UnitsStyle.DimensionLinearUnits

            Dim dimensionlinearunits30 As NXOpen.Annotations.DimensionUnit = Nothing
            dimensionlinearunits30 = sketchLinearDimensionBuilder1.Style.UnitsStyle.DimensionLinearUnits

            Dim dimensionlinearunits31 As NXOpen.Annotations.DimensionUnit = Nothing
            dimensionlinearunits31 = sketchLinearDimensionBuilder1.Style.UnitsStyle.DimensionLinearUnits

            Dim dimensionlinearunits32 As NXOpen.Annotations.DimensionUnit = Nothing
            dimensionlinearunits32 = sketchLinearDimensionBuilder1.Style.UnitsStyle.DimensionLinearUnits

            Dim dimensionlinearunits33 As NXOpen.Annotations.DimensionUnit = Nothing
            dimensionlinearunits33 = sketchLinearDimensionBuilder1.Style.UnitsStyle.DimensionLinearUnits

            Dim dimensionlinearunits34 As NXOpen.Annotations.DimensionUnit = Nothing
            dimensionlinearunits34 = sketchLinearDimensionBuilder1.Style.UnitsStyle.DimensionLinearUnits

            Dim dimensionlinearunits35 As NXOpen.Annotations.DimensionUnit = Nothing
            dimensionlinearunits35 = sketchLinearDimensionBuilder1.Style.UnitsStyle.DimensionLinearUnits

            Dim dimensionlinearunits36 As NXOpen.Annotations.DimensionUnit = Nothing
            dimensionlinearunits36 = sketchLinearDimensionBuilder1.Style.UnitsStyle.DimensionLinearUnits

            ' ----------------------------------------------
            '   Dialog Begin Linear Dimension
            ' ----------------------------------------------
            Dim markId13 As NXOpen.Session.UndoMarkId = Nothing
            markId13 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Invisible, "Linear Dimension")

            theSession.DeleteUndoMark(markId13, Nothing)

            Dim markId14 As NXOpen.Session.UndoMarkId = Nothing
            markId14 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Invisible, "Linear Dimension")

            sketchLinearDimensionBuilder1.Driving.ExpressionValue.SetFormula(r)

            sketchLinearDimensionBuilder1.Driving.ExpressionMode = NXOpen.Annotations.DrivingValueBuilder.DrivingExpressionMode.KeepExpression

            sketchLinearDimensionBuilder1.Origin.SetInferRelativeToGeometry(True)

            Dim nXObject3 As NXOpen.NXObject = Nothing
            nXObject3 = sketchLinearDimensionBuilder1.Commit()

            Dim taggedObject2 As NXOpen.TaggedObject = Nothing
            taggedObject2 = sketchLinearDimensionBuilder1.SecondAssociativity.Value

            Dim line1 As NXOpen.Line = CType(taggedObject2, NXOpen.Line)

            Dim point1_2 As NXOpen.Point3d = New NXOpen.Point3d(0.0, 0.0, 0.0)
            Dim point2_2 As NXOpen.Point3d = New NXOpen.Point3d(0.0, 0.0, 0.0)
            sketchLinearDimensionBuilder1.FirstAssociativity.SetValue(NXOpen.InferSnapType.SnapType.Start, line1, nullNXOpen_View, point1_2, Nothing, nullNXOpen_View, point2_2)

            Dim point1_3 As NXOpen.Point3d = New NXOpen.Point3d(-radius, 0.0, 0.0)
            Dim point2_3 As NXOpen.Point3d = New NXOpen.Point3d(0.0, 0.0, 0.0)
            sketchLinearDimensionBuilder1.SecondAssociativity.SetValue(NXOpen.InferSnapType.SnapType.End, line1, nullNXOpen_View, point1_3, Nothing, nullNXOpen_View, point2_3)

            sketchLinearDimensionBuilder1.Driving.ExpressionValue.SetFormula(r)

            theSession.SetUndoMarkName(markId14, "Linear Dimension - =")

            theSession.SetUndoMarkVisibility(markId14, Nothing, NXOpen.Session.MarkVisibility.Visible)

            theSession.SetUndoMarkVisibility(markId12, Nothing, NXOpen.Session.MarkVisibility.Invisible)

            Dim markId15 As NXOpen.Session.UndoMarkId = Nothing
            markId15 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Invisible, "Linear Dimension")

            Dim markId16 As NXOpen.Session.UndoMarkId = Nothing
            markId16 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Invisible, "Linear Dimension")

            Dim nXObject4 As NXOpen.NXObject = Nothing
            nXObject4 = sketchLinearDimensionBuilder1.Commit()

            theSession.DeleteUndoMark(markId16, Nothing)

            theSession.SetUndoMarkName(markId12, "Linear Dimension")

            Dim expression2 As NXOpen.Expression = sketchLinearDimensionBuilder1.Driving.ExpressionValue

            sketchLinearDimensionBuilder1.Destroy()

            theSession.DeleteUndoMark(markId15, Nothing)

            theSession.SetUndoMarkVisibility(markId12, Nothing, NXOpen.Session.MarkVisibility.Visible)

            theSession.DeleteUndoMark(markId14, Nothing)

            theSession.DeleteUndoMark(markId12, Nothing)

            Dim markId17 As NXOpen.Session.UndoMarkId = Nothing
            markId17 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Visible, "Edit Dimension")

            Dim markId18 As NXOpen.Session.UndoMarkId = Nothing
            markId18 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Invisible, "Start")

            Dim parallelDimension2 As NXOpen.Annotations.ParallelDimension = CType(theSession.ActiveSketch.FindObject("HANDLE R-3264"), NXOpen.Annotations.ParallelDimension)

            Dim sketchLinearDimensionBuilder2 As NXOpen.SketchLinearDimensionBuilder = Nothing
            sketchLinearDimensionBuilder2 = workPart.Sketches.CreateLinearDimensionBuilder(parallelDimension2)

            theSession.SetUndoMarkName(markId18, "Linear Dimension Dialog")

            theSession.SetUndoMarkVisibility(markId18, Nothing, NXOpen.Session.MarkVisibility.Visible)

            sketchLinearDimensionBuilder2.Origin.SetInferRelativeToGeometry(True)

            Dim dimensionlinearunits37 As NXOpen.Annotations.DimensionUnit = Nothing
            dimensionlinearunits37 = sketchLinearDimensionBuilder2.Style.UnitsStyle.DimensionLinearUnits

            Dim dimensionlinearunits38 As NXOpen.Annotations.DimensionUnit = Nothing
            dimensionlinearunits38 = sketchLinearDimensionBuilder2.Style.UnitsStyle.DimensionLinearUnits

            sketchLinearDimensionBuilder2.Style.OrdinateStyle.DoglegCreationOption = NXOpen.Annotations.OrdinateDoglegCreationOption.No

            Dim dimensionlinearunits39 As NXOpen.Annotations.DimensionUnit = Nothing
            dimensionlinearunits39 = sketchLinearDimensionBuilder2.Style.UnitsStyle.DimensionLinearUnits

            Dim dimensionlinearunits40 As NXOpen.Annotations.DimensionUnit = Nothing
            dimensionlinearunits40 = sketchLinearDimensionBuilder2.Style.UnitsStyle.DimensionLinearUnits

            Dim dimensionlinearunits41 As NXOpen.Annotations.DimensionUnit = Nothing
            dimensionlinearunits41 = sketchLinearDimensionBuilder2.Style.UnitsStyle.DimensionLinearUnits

            Dim dimensionlinearunits42 As NXOpen.Annotations.DimensionUnit = Nothing
            dimensionlinearunits42 = sketchLinearDimensionBuilder2.Style.UnitsStyle.DimensionLinearUnits

            sketchLinearDimensionBuilder2.Origin.SetInferRelativeToGeometry(True)

            sketchLinearDimensionBuilder2.Origin.SetInferRelativeToGeometry(True)

            sketchLinearDimensionBuilder2.Measurement.Direction = nullNXOpen_Direction

            sketchLinearDimensionBuilder2.Measurement.DirectionView = nullNXOpen_View

            sketchLinearDimensionBuilder2.Origin.SetInferRelativeToGeometry(True)

            Dim dimensionlinearunits43 As NXOpen.Annotations.DimensionUnit = Nothing
            dimensionlinearunits43 = sketchLinearDimensionBuilder2.Style.UnitsStyle.DimensionLinearUnits

            Dim dimensionlinearunits44 As NXOpen.Annotations.DimensionUnit = Nothing
            dimensionlinearunits44 = sketchLinearDimensionBuilder2.Style.UnitsStyle.DimensionLinearUnits

            Dim dimensionlinearunits45 As NXOpen.Annotations.DimensionUnit = Nothing
            dimensionlinearunits45 = sketchLinearDimensionBuilder2.Style.UnitsStyle.DimensionLinearUnits

            Dim dimensionlinearunits46 As NXOpen.Annotations.DimensionUnit = Nothing
            dimensionlinearunits46 = sketchLinearDimensionBuilder2.Style.UnitsStyle.DimensionLinearUnits

            Dim dimensionlinearunits47 As NXOpen.Annotations.DimensionUnit = Nothing
            dimensionlinearunits47 = sketchLinearDimensionBuilder2.Style.UnitsStyle.DimensionLinearUnits

            Dim dimensionlinearunits48 As NXOpen.Annotations.DimensionUnit = Nothing
            dimensionlinearunits48 = sketchLinearDimensionBuilder2.Style.UnitsStyle.DimensionLinearUnits

            Dim dimensionlinearunits49 As NXOpen.Annotations.DimensionUnit = Nothing
            dimensionlinearunits49 = sketchLinearDimensionBuilder2.Style.UnitsStyle.DimensionLinearUnits

            Dim dimensionlinearunits50 As NXOpen.Annotations.DimensionUnit = Nothing
            dimensionlinearunits50 = sketchLinearDimensionBuilder2.Style.UnitsStyle.DimensionLinearUnits

            Dim dimensionlinearunits51 As NXOpen.Annotations.DimensionUnit = Nothing
            dimensionlinearunits51 = sketchLinearDimensionBuilder2.Style.UnitsStyle.DimensionLinearUnits

            Dim dimensionlinearunits52 As NXOpen.Annotations.DimensionUnit = Nothing
            dimensionlinearunits52 = sketchLinearDimensionBuilder2.Style.UnitsStyle.DimensionLinearUnits

            Dim dimensionlinearunits53 As NXOpen.Annotations.DimensionUnit = Nothing
            dimensionlinearunits53 = sketchLinearDimensionBuilder2.Style.UnitsStyle.DimensionLinearUnits

            Dim dimensionlinearunits54 As NXOpen.Annotations.DimensionUnit = Nothing
            dimensionlinearunits54 = sketchLinearDimensionBuilder2.Style.UnitsStyle.DimensionLinearUnits

            ' ----------------------------------------------
            '   Dialog Begin Linear Dimension
            ' ----------------------------------------------
            Dim markId19 As NXOpen.Session.UndoMarkId = Nothing
            markId19 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Invisible, "Linear Dimension")

            theSession.DeleteUndoMark(markId19, Nothing)

            Dim markId20 As NXOpen.Session.UndoMarkId = Nothing
            markId20 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Invisible, "Linear Dimension")

            sketchLinearDimensionBuilder2.Driving.ExpressionValue.SetFormula(d)

            sketchLinearDimensionBuilder2.Driving.ExpressionMode = NXOpen.Annotations.DrivingValueBuilder.DrivingExpressionMode.KeepExpression

            sketchLinearDimensionBuilder2.Origin.SetInferRelativeToGeometry(True)

            Dim nXObject5 As NXOpen.NXObject = Nothing
            nXObject5 = sketchLinearDimensionBuilder2.Commit()

            Dim taggedObject3 As NXOpen.TaggedObject = Nothing
            taggedObject3 = sketchLinearDimensionBuilder2.SecondAssociativity.Value

            Dim line2 As NXOpen.Line = CType(taggedObject3, NXOpen.Line)

            Dim point1_4 As NXOpen.Point3d = New NXOpen.Point3d(0.0, 0.0, 0.0)
            Dim point2_4 As NXOpen.Point3d = New NXOpen.Point3d(0.0, 0.0, 0.0)
            sketchLinearDimensionBuilder2.FirstAssociativity.SetValue(NXOpen.InferSnapType.SnapType.Start, line2, nullNXOpen_View, point1_4, Nothing, nullNXOpen_View, point2_4)

            Dim point1_5 As NXOpen.Point3d = New NXOpen.Point3d(0.0, 0.0, -depth)
            Dim point2_5 As NXOpen.Point3d = New NXOpen.Point3d(0.0, 0.0, 0.0)
            sketchLinearDimensionBuilder2.SecondAssociativity.SetValue(NXOpen.InferSnapType.SnapType.End, line2, nullNXOpen_View, point1_5, Nothing, nullNXOpen_View, point2_5)

            sketchLinearDimensionBuilder2.Driving.ExpressionValue.SetFormula(d)

            theSession.SetUndoMarkName(markId20, "Linear Dimension - =")

            theSession.SetUndoMarkVisibility(markId20, Nothing, NXOpen.Session.MarkVisibility.Visible)

            theSession.SetUndoMarkVisibility(markId18, Nothing, NXOpen.Session.MarkVisibility.Invisible)

            Dim markId21 As NXOpen.Session.UndoMarkId = Nothing
            markId21 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Invisible, "Linear Dimension")

            Dim markId22 As NXOpen.Session.UndoMarkId = Nothing
            markId22 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Invisible, "Linear Dimension")

            Dim nXObject6 As NXOpen.NXObject = Nothing
            nXObject6 = sketchLinearDimensionBuilder2.Commit()

            theSession.DeleteUndoMark(markId22, Nothing)

            theSession.SetUndoMarkName(markId18, "Linear Dimension")

            Dim expression3 As NXOpen.Expression = sketchLinearDimensionBuilder2.Driving.ExpressionValue

            sketchLinearDimensionBuilder2.Destroy()

            theSession.DeleteUndoMark(markId21, Nothing)

            theSession.SetUndoMarkVisibility(markId18, Nothing, NXOpen.Session.MarkVisibility.Visible)

            theSession.DeleteUndoMark(markId20, Nothing)

            theSession.DeleteUndoMark(markId18, Nothing)

            Dim scaleAboutPoint8 As NXOpen.Point3d = New NXOpen.Point3d(87.409432970314015, -1.2955943622082449, 0.0)
            Dim viewCenter8 As NXOpen.Point3d = New NXOpen.Point3d(-87.409432970313901, 1.295594362208156, 0.0)
            workPart.ModelingViews.WorkView.ZoomAboutPoint(1.25, scaleAboutPoint8, viewCenter8)

            ' ----------------------------------------------
            '   Menu: Task->Finish Sketch
            ' ----------------------------------------------
            Dim revolve1 As NXOpen.Features.Revolve = CType(workPart.Features.FindObject("REVOLVED(2)"), NXOpen.Features.Revolve)

            Dim section1() As NXOpen.Section
            section1 = revolve1.GetSections()

            section1(0).CleanMappingData()

            theSession.Preferences.Sketch.SectionView = False

            Dim markId23 As NXOpen.Session.UndoMarkId = Nothing
            markId23 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Invisible, "Deactivate Sketch")

            theSession.ActiveSketch.Deactivate(NXOpen.Sketch.ViewReorient.True, NXOpen.Sketch.UpdateLevel.Model)

            theSession.DeleteUndoMarksSetInTaskEnvironment()

            theSession.EndTaskEnvironment()

            ' ----------------------------------------------
            '   Menu: File->Export->STL...
            ' ----------------------------------------------
            Dim markId24 As NXOpen.Session.UndoMarkId = Nothing
            markId24 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Visible, "Start")

            Dim sTLCreator1 As NXOpen.STLCreator = Nothing
            sTLCreator1 = theSession.DexManager.CreateStlCreator()

            sTLCreator1.AutoNormalGen = True

            sTLCreator1.ChordalTol = 0.080000000000000002

            sTLCreator1.AdjacencyTol = 0.080000000000000002

            theSession.SetUndoMarkName(markId24, "STL Export Dialog")

            Dim body1 As NXOpen.Body = CType(workPart.Bodies.FindObject("REVOLVED(2)"), NXOpen.Body)

            Dim added1 As Boolean = Nothing
            added1 = sTLCreator1.ExportSelectionBlock.Add(body1)

            Dim markId25 As NXOpen.Session.UndoMarkId = Nothing
            markId25 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Invisible, "STL Export")

            theSession.DeleteUndoMark(markId25, Nothing)

            Dim markId26 As NXOpen.Session.UndoMarkId = Nothing
            markId26 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Invisible, "STL Export")

            sTLCreator1.OutputFile = "\\ad.uillinois.edu\engr-ews\ibansal2\Desktop\Plumes_Research\Test Crater Type 2 STLs\r" + r + "_d" + d + "_c" + c + "_t2.stl"

            Dim nXObject7 As NXOpen.NXObject = Nothing
            nXObject7 = sTLCreator1.Commit()

            theSession.DeleteUndoMark(markId26, Nothing)

            theSession.SetUndoMarkName(markId24, "STL Export")

            sTLCreator1.Destroy()
        Next
    Next
Next

' ----------------------------------------------
'   Menu: Tools->Journal->Stop Recording
' ----------------------------------------------

End Sub
End Module
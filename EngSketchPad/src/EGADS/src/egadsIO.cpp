/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             Load & Save Functions
 *
 *      Copyright 2011-2025, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

//#define WRITECSYS

#include "egadsTypes.h"
#include "egadsInternals.h"
#include "egadsClasses.h"
#include <IGESControl_Controller.hxx>
#include <IGESData_IGESModel.hxx>
#include <IGESData_Protocol.hxx>
#include <IGESBasic_Name.hxx>
#include <IGESGraph_Color.hxx>
#include <IGESSelect_WorkLibrary.hxx>
#include <IGESCAFControl.hxx>
#include <IGESSolid_Face.hxx>
#include <IGESSolid_EdgeList.hxx>
#include <IGESGeom_CompositeCurve.hxx>
#include <TransferBRep_ShapeMapper.hxx>
#ifdef STEPASSATTRS
#include <StepBasic_Product.hxx>
#include <StepBasic_ProductDefinition.hxx>
#include <StepBasic_ProductDefinitionFormation.hxx>
#include <StepRepr_ProductDefinitionShape.hxx>
#include <StepRepr_NextAssemblyUsageOccurrence.hxx>
#endif
#include <StepRepr_RepresentationItem.hxx>
#include <StepRepr_RepresentationContext.hxx>
#include <StepRepr_ShapeRepresentationRelationship.hxx>
#ifdef WRITECSYS
#include <StepRepr_ConstructiveGeometryRepresentation.hxx>
#include <StepRepr_HArray1OfRepresentationItem.hxx>
#include <StepGeom_Axis2Placement3d.hxx>
#include <GeomToStep_MakeAxis2Placement3d.hxx>
#endif
#include <STEPConstruct.hxx>
#include <STEPConstruct_Styles.hxx>
#include <STEPControl_ActorWrite.hxx>
#include <StepData_StepModel.hxx>
#include <StepShape_VertexPoint.hxx>
#include <StepShape_EdgeCurve.hxx>
#include <StepShape_FaceSurface.hxx>
#include <StepGeom_Point.hxx>
#include <StepGeom_Curve.hxx>
#include <StepGeom_TrimmedCurve.hxx>
#include <StepGeom_Surface.hxx>
#include <StepGeom_SurfaceCurve.hxx>
#include <StepVisual_Colour.hxx>
#include <StepVisual_ColourRgb.hxx>
#include <StepVisual_PreDefinedColour.hxx>
#include <StepVisual_PreDefinedItem.hxx>
#include <StepVisual_DraughtingPreDefinedColour.hxx>
#include <StepVisual_StyledItem.hxx>
#include <StepVisual_StyledItemTarget.hxx>
#include <StepVisual_OverRidingStyledItem.hxx>
#include <StepVisual_PresentationRepresentation.hxx>
#include <StepVisual_ContextDependentOverRidingStyledItem.hxx>
#include <StepVisual_MechanicalDesignGeometricPresentationRepresentation.hxx>
#include <StepVisual_PresentationLayerAssignment.hxx>
#include <StepShape_Shell.hxx>
#include <StepShape_OpenShell.hxx>
#include <StepShape_ClosedShell.hxx>
#include <StepShape_ShapeRepresentation.hxx>
#include <StepShape_ShellBasedSurfaceModel.hxx>
#include <StepShape_ContextDependentShapeRepresentation.hxx>
#include <StepShape_GeometricCurveSet.hxx>
#include <StepShape_GeometricSetSelect.hxx>
#include <Transfer_TransientListBinder.hxx>
#include <Quantity_Color.hxx>
#include <Quantity_ColorRGBA.hxx>
#include <Quantity_NameOfColor.hxx>
#include <Interface_Graph.hxx>
#include <Interface_EntityIterator.hxx>
#include <XSControl_WorkSession.hxx>
#include <XSControl_TransferReader.hxx>
#include <XSControl_Utils.hxx>
#include <Transfer_ResultFromModel.hxx>
#include <Transfer_TransientProcess.hxx>
#include <Transfer_ResultFromTransient.hxx>
#include <TransferBRep.hxx>
#include <TransferBRep_ShapeBinder.hxx>
#include <TransferBRep_ShapeListBinder.hxx>
#include <MoniTool_Macros.hxx>
#include <MoniTool_DataMapOfShapeTransient.hxx>
#include <APIHeaderSection_MakeHeader.hxx>
#include <Interface_Static.hxx>


#ifdef WIN32
#define snprintf _snprintf
#endif


//#define INTERIM

#define UVTOL    1.e-4


  extern "C" void EG_revision( int *major, int *minor, char **OCCrev );
  extern "C" void EG_initOCC( );
  extern "C" int  EG_destroyTopology( egObject *topo );
  extern "C" int  EG_fullAttrs( const egObject *obj );
  extern "C" void EG_attrBuildSeq( egAttrs *attrs );
  extern "C" void EG_readAttrs( egObject *obj, int nattr, FILE *fp );
  extern "C" void EG_writeAttr( egAttrs *attrs, FILE *fp );
  extern "C" int  EG_writeNumAttr( egAttrs *attrs );
  extern "C" int  EG_dereferenceTopObj( egObject *object,
                                        /*@null@*/ const egObject *ref );

  extern "C" int  EG_loadModel( egObject *context, int bflg, const char *name,
                                egObject **model );
  extern "C" int  EG_saveModel( const egObject *model, const char *name );
  extern "C" int  EG_saveTess( egObject *tess, const char *name );
  extern "C" int  EG_loadTess( egObject *body, const char *name,
                               egObject **tess );

  extern "C" int  EG_attributeAdd( egObject *obj, const char *name, int type,
                                   int len, /*@null@*/ const int    *ints,
                                            /*@null@*/ const double *reals,
                                            /*@null@*/ const char   *str );
  extern "C" int  EG_attributeRet( const egObject *obj, const char *name,
                                   int *type, int *len,
                                   /*@null@*/ const int    **ints,
                                   /*@null@*/ const double **reals,
                                   /*@null@*/ const char   **str );
#ifdef WRITECSYS
  extern "C" int  EG_attributeGet( const egObject *obj, int index,
                                   const char **name, int *type, int *len,
                                   /*@null@*/ const int    **ints,
                                   /*@null@*/ const double **reals,
                                   /*@null@*/ const char   **str );
  extern "C" int  EG_attributeNum( const egObject *obj, int *num );
#endif
  extern "C" int  EG_statusTessBody( egObject *tess, egObject **body,
                                     int *state, int *npts );
  extern "C" int  EG_getBodyTopos( const egObject *body, egObject *src,
                                   int oclass, int *ntopo, egObject ***topos );
  extern "C" int  EG_getTessEdge( const egObject *tess, int indx, int *len,
                                  const double **xyz, const double **t );
  extern "C" int  EG_getTessFace( const egObject *tess, int indx, int *len,
                                  const double **xyz, const double **uv,
                                  const int **ptype, const int **pindex,
                                  int *ntri, const int **tris,
                                  const int **tric );
  extern "C" int  EG_initTessBody( egObject *object, egObject **tess );
  extern "C" int  EG_setTessEdge( const egObject *tess, int index, int len,
                                  const double *xyz, const double *t );
  extern "C" int  EG_setTessFace( const egObject *tess, int index, int len,
                                  const double *xyz, const double *uv,
                                  int ntri, const int *tris );
  extern "C" int  EG_isSame ( const egObject *obj1, const egObject *obj2 );
  extern "C" int  EG_invEvaluate( const egObject *obj, double *xyz,
                                  double *param, double *results );
  extern "C" int  EG_writeEBody( const egObject *EBody, FILE *fp );
  extern "C" int  EG_readEBody( FILE *fp, egObject *body, egObject **EBody );

  extern     void EG_splitPeriodics( int ib, egadsBody *body,
                                     egadsShapeData &labels );
  extern     void EG_splitMultiplicity( egadsBody *body, egadsShapeData &labels,
                                        int outLevel );
  extern     int  EG_traverseBody( egObject *context, int i, egObject *bobj,
                                   egObject *topObj, egadsBody *body,
                                   int *nerr );


void
EG_revision(int *major, int *minor, char **OCCrev)
{
#ifdef INTERIM
  static char OCCrevStr[42];

  snprintf(OCCrevStr, 41, "Interim Release with OpenCASCADE %d.%d.%d",
           OCC_VERSION_MAJOR, OCC_VERSION_MINOR, OCC_VERSION_MAINTENANCE);
#else
  static char OCCrevStr[26];

  snprintf(OCCrevStr, 25, "with OpenCASCADE %d.%d.%d", OCC_VERSION_MAJOR,
           OCC_VERSION_MINOR, OCC_VERSION_MAINTENANCE);
#endif
  *major  = EGADSMAJOR;
  *minor  = EGADSMINOR;
  *OCCrev = OCCrevStr;
}


TopoDS_Shape
egadsShapeData::Update(const TopoDS_Shape& oldShape,
                       const TopoDS_Shape& newShape)
{
  Standard_Integer i;
  if (labels.Extent() == 0 &&
      colors.Extent() == 0 &&
      bodycolors.Extent() == 0) return newShape;

  i = labels.FindIndex(oldShape);
  if (i > 0)
    labels.Add(newShape, labels(i));

  i = colors.FindIndex(oldShape);
  if (i > 0)
    colors.Add(newShape, colors(i));

  i = bodycolors.FindIndex(oldShape);
  if (i > 0)
    bodycolors.Add(newShape, bodycolors(i));

  TopTools_IndexedMapOfShape oldMap;
  TopExp::MapShapes(oldShape, oldMap);

  TopTools_IndexedMapOfShape newMap;
  TopExp::MapShapes(newShape, newMap);

  if (oldMap.Extent() != newMap.Extent()) return newShape;

  for (int j = 1; j <= oldMap.Extent(); j++) {
    TopoDS_Shape shape = oldMap(j);
    if (shape == newMap(j)) continue;
    i = labels.FindIndex(shape);
    if (i > 0)
      labels.Add(newMap(j), labels(i));
    i = colors.FindIndex(shape);
    if (i > 0)
      colors.Add(newMap(j), colors(i));
    i = bodycolors.FindIndex(shape);
    if (i > 0)
      bodycolors.Add(newMap(j), bodycolors(i));
  }

  return newShape;
}


TopoDS_Shape
egadsShapeData::Update(const TopoDS_Shape& oldShape,
                       BRepBuilderAPI_ModifyShape& xForm)
{
  TopoDS_Shape newShape = xForm.ModifiedShape(oldShape);

  Standard_Integer i;
  if (labels.Extent() == 0 &&
      colors.Extent() == 0 &&
      bodycolors.Extent() == 0) return newShape;

  i = labels.FindIndex(oldShape);
  if (i > 0)
    labels.Add(newShape, labels(i));

  i = colors.FindIndex(oldShape);
  if (i > 0)
    colors.Add(newShape, colors(i));

  i = bodycolors.FindIndex(oldShape);
  if (i > 0)
    bodycolors.Add(newShape, bodycolors(i));

  TopTools_IndexedMapOfShape oldMap;
  TopExp::MapShapes(oldShape, oldMap);

  for (int j = 1; j <= oldMap.Extent(); j++) {

    TopoDS_Shape shape = oldMap(j);

    const TopTools_ListOfShape& mods = xForm.Modified(shape);

    if (mods.Extent() == 0) continue;
    i = labels.FindIndex(shape);
    if (i > 0) {
      const Label& label = labels(i);

      TopTools_ListIteratorOfListOfShape it(mods);
      for (; it.More(); it.Next()) {
        labels.Add(it.Value(), label);
      }
    }

    i = colors.FindIndex(shape);
    if (i > 0) {
      const Quantity_Color& color = colors(i);

      TopTools_ListIteratorOfListOfShape it(mods);
      for (; it.More(); it.Next()) {
        colors.Add(it.Value(), color);
      }
    }

    i = bodycolors.FindIndex(shape);
    if (i > 0) {
      const BodyColor& bodycolor = bodycolors(i);

      TopTools_ListIteratorOfListOfShape it(mods);
      for (; it.More(); it.Next()) {
        bodycolors.Add(it.Value(), bodycolor);
      }
    }
  }

  return newShape;
}


TopoDS_Shape
egadsShapeData::Update(const TopoDS_Shape& oldShape,
                       const Handle(BRepTools_ReShape)& reShape)
{
  TopoDS_Shape newShape = reShape->Apply(oldShape);

  return Update(oldShape, newShape, reShape->History());
}


TopoDS_Shape
egadsShapeData::Update(const TopoDS_Shape& oldShape,
                       const TopoDS_Shape& newShape,
                       const Handle(BRepTools_History)& history)
{
  Standard_Integer i;
  if (labels.Extent() == 0 && 
      colors.Extent() == 0 &&
      bodycolors.Extent() == 0) return newShape;

  i = labels.FindIndex(oldShape);
  if (i > 0)
    labels.Add(newShape, labels(i));

  i = colors.FindIndex(oldShape);
  if (i > 0)
    colors.Add(newShape, colors(i));

  i = bodycolors.FindIndex(oldShape);
  if (i > 0)
    bodycolors.Add(newShape, bodycolors(i));

  TopTools_IndexedMapOfShape oldMap;
  TopExp::MapShapes(oldShape, oldMap);

  for (int j = 1; j <= oldMap.Extent(); j++) {

    TopoDS_Shape shape = oldMap(j);

    const TopTools_ListOfShape& mods = history->Modified(shape);

    if (mods.Extent() > 0) {
      i = labels.FindIndex(shape);
      if (i > 0) {
        const Label& label = labels(i);

        TopTools_ListIteratorOfListOfShape it(mods);
        for (; it.More(); it.Next()) {
          labels.Add(it.Value(), label);
        }
      }

      i = colors.FindIndex(shape);
      if (i > 0) {
        const Quantity_Color& color = colors(i);

        TopTools_ListIteratorOfListOfShape it(mods);
        for (; it.More(); it.Next()) {
          colors.Add(it.Value(), color);
        }
      }

      i = bodycolors.FindIndex(shape);
      if (i > 0) {
        const BodyColor& bodycolor = bodycolors(i);

        TopTools_ListIteratorOfListOfShape it(mods);
        for (; it.More(); it.Next()) {
          bodycolors.Add(it.Value(), bodycolor);
        }
      }
    }

    const TopTools_ListOfShape& gens = history->Generated(shape);

    if (gens.Extent() > 0) {
      i = labels.FindIndex(shape);
      if (i > 0) {
        const Label& label = labels(i);

        TopTools_ListIteratorOfListOfShape it(gens);
        for (; it.More(); it.Next()) {
          labels.Add(it.Value(), label);
        }
      }

      i = colors.FindIndex(shape);
      if (i > 0) {
        const Quantity_Color& color = colors(i);

        TopTools_ListIteratorOfListOfShape it(gens);
        for (; it.More(); it.Next()) {
          colors.Add(it.Value(), color);
        }
      }

      i = bodycolors.FindIndex(shape);
      if (i > 0) {
        const BodyColor& bodycolor = bodycolors(i);

        TopTools_ListIteratorOfListOfShape it(gens);
        for (; it.More(); it.Next()) {
          bodycolors.Add(it.Value(), bodycolor);
        }
      }
    }
  }

  return newShape;
}


void
egadsShapeData::AddShapeName(egadsMap& shapes)
{
  int index;
  for (int i = 1; i <= shapes.map.Extent(); i++) {
    const TopoDS_Shape& shape = shapes.map.FindKey(i);
    
    index = labels.FindIndex(shape);
    if (index > 0) {
      const char *name = labels(index).shapeName;
      EG_attributeAdd(shapes.objs[i-1], "Name", ATTRSTRING, 1, NULL, NULL, name);
    }
  }
}


void
egadsShapeData::AddShapeColor(egadsMap& shapes)
{
  int index;
  for (int i = 1; i <= shapes.map.Extent(); i++) {
    const TopoDS_Shape& shape = shapes.map.FindKey(i);

    index = colors.FindIndex(shape);
    if (index > 0) {
      const Quantity_Color& color = colors(index);

      const double rgb[3] = {color.Red(), color.Green(), color.Blue()};

      const char* name = NULL;
      if (rgb[0] == 1   && rgb[1] == 0   && rgb[2] == 0  ) name = "red";
      else if (rgb[0] == 0   && rgb[1] == 1   && rgb[2] == 0  ) name = "green";
      else if (rgb[0] == 0   && rgb[1] == 0   && rgb[2] == 1  ) name = "blue";
      else if (rgb[0] == 1   && rgb[1] == 1   && rgb[2] == 0  ) name = "yellow";
      else if (rgb[0] == 1   && rgb[1] == 0   && rgb[2] == 1  ) name = "magenta";
      else if (rgb[0] == 0   && rgb[1] == 1   && rgb[2] == 1  ) name = "cyan";
      else if (rgb[0] == 1   && rgb[1] == 1   && rgb[2] == 1  ) name = "white";
      else if (rgb[0] == 0   && rgb[1] == 0   && rgb[2] == 0  ) name = "black";
      else if (rgb[0] == 1   && rgb[1] == 0.5 && rgb[2] == 0.5) name = "lred";
      else if (rgb[0] == 0.5 && rgb[1] == 1   && rgb[2] == 0.5) name = "lgreen";
      else if (rgb[0] == 0.5 && rgb[1] == 0.5 && rgb[2] == 1  ) name = "lblue";

      if (name != NULL) {
        EG_attributeAdd(shapes.objs[i-1],  "Color", ATTRSTRING, 1, NULL, NULL,
                        name);
      } else {
        EG_attributeAdd(shapes.objs[i-1],  "Color", ATTRREAL  , 3, NULL, rgb ,
                        NULL);
      }
    }
  }
}


void
EG_attriBodyTrav(const egObject *obj, const TopoDS_Shape& shape,
                 egadsBody *pbody)
{
  if (obj->blind == NULL) return;

  if (obj->oclass == NODE) {

    int index = pbody->nodes.map.FindIndex(shape);
    if (index != 0)
      EG_attributeDup(obj, pbody->nodes.objs[index-1]);
    else
      printf(" EGADS Internal: Dropping Node attributes!\n");

  } else if (obj->oclass == EDGE) {

    egadsEdge *pedge = (egadsEdge *) obj->blind;
    int index = pbody->edges.map.FindIndex(shape);
    if (index != 0) {
      EG_attributeDup(obj, pbody->edges.objs[index-1]);
    } else {
      printf(" EGADS Internal: Dropping Edge attributes!\n");
    }

    TopoDS_Vertex V1, V2;
    TopExp::Vertices(TopoDS::Edge(shape), V1, V2);

    EG_attriBodyTrav(pedge->nodes[0], V1, pbody);
    if (obj->mtype == TWONODE)
      EG_attriBodyTrav(pedge->nodes[1], V2, pbody);

  } else if (obj->oclass == LOOP) {

    egadsLoop *ploop = (egadsLoop *) obj->blind;
    int index = pbody->loops.map.FindIndex(shape);
    if (index != 0) {
      EG_attributeDup(obj, pbody->loops.objs[index-1]);
    } else {
      printf(" EGADS Internal: Dropping Loop attributes!\n");
    }

    BRepTools_WireExplorer Exp(ploop->loop);
    for (int i = 0; Exp.More(); Exp.Next(), i++) {
      EG_attriBodyTrav(ploop->edges[i], Exp.Current(), pbody);
    }

  } else if (obj->oclass == FACE) {

    egadsFace *pface = (egadsFace *) obj->blind;
    int index = pbody->faces.map.FindIndex(shape);
    if (index != 0) {
      EG_attributeDup(obj, pbody->faces.objs[index-1]);
    } else {
      printf(" EGADS Internal: Dropping Face attributes!\n");
    }

    TopExp_Explorer Exp(pface->face, TopAbs_WIRE);
    for (int i = 0; Exp.More(); Exp.Next(), i++) {
      EG_attriBodyTrav(pface->loops[i], Exp.Current(), pbody);
    }

  } else if (obj->oclass == SHELL) {

    egadsShell *pshell = (egadsShell *) obj->blind;
    int index = pbody->shells.map.FindIndex(shape);
    if (index != 0) {
      EG_attributeDup(obj, pbody->shells.objs[index-1]);
    } else {
      printf(" EGADS Internal: Dropping Shell attributes!\n");
    }

    TopExp_Explorer Exp(pshell->shell, TopAbs_FACE);
    for (int i = 0; Exp.More(); Exp.Next(), i++) {
      EG_attriBodyTrav(pshell->faces[i], Exp.Current(), pbody);
    }
  }
}


int
EG_attriBodyDup(const egObject *src, egObject *dst)
{
  int          i, j, n, stat, nents, nattr;
  double       tmin, tmax, trange[2], xyz[3];
  egObject     *aobj, *dobj, *geom;
  egAttrs      *attrs;
  TopoDS_Shape shape;

  if ((src == NULL) || (dst == NULL)) return EGADS_NULLOBJ;
  if  (src->magicnumber != MAGIC)     return EGADS_NOTOBJ;
  if  (src->oclass < NODE)            return EGADS_NOTTOPO;
  if  (src->blind == NULL)            return EGADS_NODATA;
  int outLevel = EG_outLevel(src);
  int fullAttr = EG_fullAttrs(src);

  if (src->oclass == MODEL) {
    egadsModel *pmdl = (egadsModel *) src->blind;
    for (i = 0; i < pmdl->nbody; i++) {
      stat = EG_attriBodyDup(pmdl->bodies[i], dst);
      if (stat != EGADS_SUCCESS) return stat;
    }
    return EGADS_SUCCESS;
  }
  if (dst->magicnumber != MAGIC) {
    if (outLevel > 0)
      printf(" EGADS Error: dst not an EGO (EG_attriBodyDup)!\n");
    return EGADS_NOTOBJ;
  }
  if (dst->oclass != BODY) {
    if (outLevel > 0)
      printf(" EGADS Error: dst not a BODY (EG_attriBodyDup)!\n");
    return EGADS_NOTBODY;
  }
  if (dst->blind == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL dst pointer (EG_attriBodyDup)!\n");
    return EGADS_NODATA;
  }
  egadsBody *pbody = (egadsBody *) dst->blind;

  if (src->oclass == BODY) {

    // use hashed body data on both ends
    egadsBody *pbods = (egadsBody *) src->blind;
    shape = pbody->shape;
    if (shape.IsSame(pbods->shape)) EG_attributeDup(src, dst);

    nents = pbods->shells.map.Extent();
    for (i = 0; i < nents; i++) {
      aobj = pbods->shells.objs[i];
      if (aobj->attrs == NULL) continue;
      attrs = (egAttrs *) aobj->attrs;
      nattr = attrs->nattrs;
      if (nattr <= 0) continue;
      shape = pbods->shells.map(i+1);
      j     = pbody->shells.map.FindIndex(shape);
      if (j == 0) continue;                     // not in the dst body
      dobj = pbody->shells.objs[j-1];
      EG_attributeDup(aobj, dobj);
    }

    nents = pbods->faces.map.Extent();
    for (i = 0; i < nents; i++) {
      aobj = pbods->faces.objs[i];
      if (aobj->attrs == NULL) continue;
      attrs = (egAttrs *) aobj->attrs;
      nattr = attrs->nattrs;
      if (nattr <= 0) continue;
      shape = pbods->faces.map(i+1);
      j     = pbody->faces.map.FindIndex(shape);
      if (j == 0) continue;                     // not in the dst body
      dobj = pbody->faces.objs[j-1];
      EG_attributeDup(aobj, dobj);
    }

    nents = pbods->loops.map.Extent();
    for (i = 0; i < nents; i++) {
      aobj = pbods->loops.objs[i];
      if (aobj->attrs == NULL) continue;
      attrs = (egAttrs *) aobj->attrs;
      nattr = attrs->nattrs;
      if (nattr <= 0) continue;
      shape = pbods->loops.map(i+1);
      j     = pbody->loops.map.FindIndex(shape);
      if (j == 0) continue;                     // not in the dst body
      dobj = pbody->loops.objs[j-1];
      EG_attributeDup(aobj, dobj);
    }

    nents = pbods->edges.map.Extent();
    for (i = 0; i < nents; i++) {
      aobj = pbods->edges.objs[i];
      if (aobj->attrs == NULL) continue;
      attrs = (egAttrs *) aobj->attrs;
      nattr = attrs->nattrs;
      if (nattr <= 0) continue;
      shape = pbods->edges.map(i+1);
      j     = pbody->edges.map.FindIndex(shape);
      if (j == 0) {
        if (fullAttr == 0) continue;
        if (aobj->mtype == DEGENERATE) continue;
        egadsEdge *pedge = (egadsEdge *) aobj->blind;
        trange[0] = pedge->trange[0];
        trange[1] = pedge->trange[1];
        geom      = pedge->curve;
        n = pbody->edges.map.Extent();
        for (j = 1; j <= n; j++) {
          dobj = pbody->edges.objs[j-1];
          if (dobj->mtype == DEGENERATE) continue;
          if (EG_isSame(aobj, dobj) == 0) {
            if (dobj->mtype == ONENODE) {
              if (aobj->mtype == ONENODE) EG_attributeDup(aobj, dobj);
              continue;
            }
            egadsEdge *pedgd = (egadsEdge *) dobj->blind;
            egadsNode *pnod0 = (egadsNode *) pedgd->nodes[0]->blind;
            egadsNode *pnod1 = (egadsNode *) pedgd->nodes[1]->blind;
            stat = EG_invEvaluate(geom, pnod0->xyz, &tmin, xyz);
            if (stat != EGADS_SUCCESS) continue;
            stat = EG_invEvaluate(geom, pnod1->xyz, &tmax, xyz);
            if (stat != EGADS_SUCCESS) continue;
            if ((tmin+UVTOL < trange[0]) || (tmax-UVTOL > trange[1])) continue;
            EG_attributeDup(aobj, dobj);
          }
        }
      } else {
        dobj = pbody->edges.objs[j-1];
        EG_attributeDup(aobj, dobj);
      }
    }

    nents = pbods->nodes.map.Extent();
    for (i = 0; i < nents; i++) {
      aobj = pbods->nodes.objs[i];
      if (aobj->attrs == NULL) continue;
      attrs = (egAttrs *) aobj->attrs;
      nattr = attrs->nattrs;
      if (nattr <= 0) continue;
      shape = pbods->nodes.map(i+1);
      j     = pbody->nodes.map.FindIndex(shape);
      if (j == 0) {
        if (fullAttr == 0) continue;
        n = pbody->nodes.map.Extent();
        for (j = 1; j <= n; j++) {
          dobj = pbody->nodes.objs[j-1];
          if (EG_isSame(aobj, dobj) == 0) {
            EG_attributeDup(aobj, dobj);
            break;
          }
        }
      } else {
        dobj = pbody->nodes.objs[j-1];
        EG_attributeDup(aobj, dobj);
      }
    }

  } else {

    // traverse the source to find objects with attributes
    if (((dst->mtype == SOLIDBODY) || (dst->mtype == SHEETBODY)) &&
         (src->oclass == SHELL)) {
      egadsShell *pshell = (egadsShell *) src->blind;
      TopExp_Explorer Exp(pbody->shape, TopAbs_SHELL);
      bool found = false;
      for (; Exp.More(); Exp.Next()) {
        if (pshell->shell.IsSame(Exp.Current())) {
          EG_attriBodyTrav(src, Exp.Current(), pbody);
          found = true;
          break;
        }
      }
      if (!found)
        printf(" EGADS Internal: Dropping Shell and sub-shape attributes!\n");
    } else {
      EG_attriBodyTrav(src, pbody->shape, pbody);
    }
  }

  return EGADS_SUCCESS;
}


int
EG_attriBodyCopy(const egObject *src, /*@null@*/ double *xform, egObject *dst)
{
  int      i, nents, nattr;
  egObject *aobj, *dobj;
  egAttrs  *attrs;

  if ((src == NULL) || (dst == NULL)) return EGADS_NULLOBJ;
  if  (src->magicnumber != MAGIC)     return EGADS_NOTOBJ;
  if  (src->oclass < NODE)            return EGADS_NOTTOPO;
  if  (src->blind == NULL)            return EGADS_NODATA;
  int outLevel = EG_outLevel(src);

  if (src->oclass == MODEL) {
    if (outLevel > 0)
      printf(" EGADS Error: src MODEL not supported (EG_attriBodyCopy)!\n");
    return EGADS_NOTMODEL;
  }
  if (dst->magicnumber != MAGIC) {
    if (outLevel > 0)
      printf(" EGADS Error: dst not an EGO (EG_attriBodyCopy)!\n");
    return EGADS_NOTOBJ;
  }
  if (src->oclass != BODY) {
    if (outLevel > 0)
      printf(" EGADS Error: src not a BODY (EG_attriBodyCopy)!\n");
    return EGADS_NOTBODY;
  }
  if (dst->oclass != BODY) {
    if (outLevel > 0)
      printf(" EGADS Error: dst not a BODY (EG_attriBodyCopy)!\n");
    return EGADS_NOTBODY;
  }
  if (src->blind == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL src pointer (EG_attriBodyCopy)!\n");
    return EGADS_NODATA;
  }
  if (dst->blind == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL dst pointer (EG_attriBodyCopy)!\n");
    return EGADS_NODATA;
  }
  EG_attributeXDup(src, xform, dst);
  egadsBody *pbods = (egadsBody *) src->blind;
  egadsBody *pbody = (egadsBody *) dst->blind;

  nents = pbods->shells.map.Extent();
  for (i = 0; i < nents; i++) {
    aobj = pbods->shells.objs[i];
    if (aobj->attrs == NULL) continue;
    attrs = (egAttrs *) aobj->attrs;
    nattr = attrs->nattrs;
    if (nattr <= 0) continue;
    dobj = pbody->shells.objs[i];
    EG_attributeXDup(aobj, xform, dobj);
  }

  nents = pbods->faces.map.Extent();
  for (i = 0; i < nents; i++) {
    aobj = pbods->faces.objs[i];
    if (aobj->attrs == NULL) continue;
    attrs = (egAttrs *) aobj->attrs;
    nattr = attrs->nattrs;
    if (nattr <= 0) continue;
    dobj = pbody->faces.objs[i];
    EG_attributeXDup(aobj, xform, dobj);
  }

  nents = pbods->loops.map.Extent();
  for (i = 0; i < nents; i++) {
    aobj = pbods->loops.objs[i];
    if (aobj->attrs == NULL) continue;
    attrs = (egAttrs *) aobj->attrs;
    nattr = attrs->nattrs;
    if (nattr <= 0) continue;
    dobj = pbody->loops.objs[i];
    EG_attributeXDup(aobj, xform, dobj);
  }

  nents = pbods->edges.map.Extent();
  for (i = 0; i < nents; i++) {
    aobj = pbods->edges.objs[i];
    if (aobj->attrs == NULL) continue;
    attrs = (egAttrs *) aobj->attrs;
    nattr = attrs->nattrs;
    if (nattr <= 0) continue;
    dobj = pbody->edges.objs[i];
    EG_attributeXDup(aobj, xform, dobj);
  }

  nents = pbods->nodes.map.Extent();
  for (i = 0; i < nents; i++) {
    aobj = pbods->nodes.objs[i];
    if (aobj->attrs == NULL) continue;
    attrs = (egAttrs *) aobj->attrs;
    nattr = attrs->nattrs;
    if (nattr <= 0) continue;
    dobj = pbody->nodes.objs[i];
    EG_attributeXDup(aobj, xform, dobj);
  }

  return EGADS_SUCCESS;
}


void
EG_readAttrs(egObject *obj, int nattr, FILE *fp)
{
  int     i, j, n, type, namlen, len, ival, nseq, *ivec = NULL;
  char    *name, cval, *string = NULL;
  double  rval, *rvec = NULL;
  egAttrs *attrs = NULL;
  egAttr  *attr  = NULL;

  attr = (egAttr *) EG_alloc(nattr*sizeof(egAttr));
  if (attr != NULL) {
    attrs = (egAttrs *) EG_alloc(sizeof(egAttrs));
    if (attrs == NULL) {
      EG_free(attr);
      attr = NULL;
    }
  }

  for (nseq = n = i = 0; i < nattr; i++) {
    j = fscanf(fp, "%d %d %d", &type, &namlen, &len);
    if (j != 3) break;
    name = NULL;
    if ((attrs != NULL) && (namlen != 0))
      name = (char *) EG_alloc((namlen+1)*sizeof(char));
    if (name != NULL) {
      fscanf(fp, "%s", name);
      for (j = 0; j < namlen; j++)
        if (name[j] == 127) {
          nseq++;
          name[j] = 32;
          break;
        }
    } else {
      for (j = 0; j < namlen; j++) fscanf(fp, "%c", &cval);
    }
    if (type == ATTRINT) {
      if (len == 1) {
        fscanf(fp, "%d", &ival);
      } else {
        ivec = NULL;
        if ((name != NULL) && (len != 0))
          ivec = (int *) EG_alloc(len*sizeof(int));
        if (ivec == NULL) {
          for (j = 0; j < len; j++) fscanf(fp, "%d", &ival);
          if (name != NULL) {
            EG_free(name);
            name = NULL;
          }
        } else {
          for (j = 0; j < len; j++) fscanf(fp, "%d", &ivec[j]);
        }
      }
    } else if ((type == ATTRREAL) || (type == ATTRCSYS)) {
      if (len == 1) {
        fscanf(fp, "%lf", &rval);
      } else {
        rvec = NULL;
        if ((name != NULL) && (len != 0))
          rvec = (double *) EG_alloc(len*sizeof(double));
        if (rvec == NULL) {
          for (j = 0; j < len; j++) fscanf(fp, "%lf", &rval);
          if (name != NULL) {
            EG_free(name);
            name = NULL;
          }
        } else {
          for (j = 0; j < len; j++) fscanf(fp, "%lf", &rvec[j]);
        }
      }
    } else {
      do {
        j = fscanf(fp, "%c", &cval);
        if (j == 0) break;
      } while (cval != '#');
      string = NULL;
      if (name != NULL)
        string = (char *) EG_alloc((len+1)*sizeof(char));
      if (string != NULL) {
        for (j = 0; j < len; j++) fscanf(fp, "%c", &string[j]);
        string[len] = 0;
      } else {
        for (j = 0; j < len; j++) fscanf(fp, "%c", &cval);
        if (name != NULL) {
          EG_free(name);
          name = NULL;
        }
      }
    }

    if (name != NULL) {
      attr[n].name   = name;
      attr[n].type   = type;
      attr[n].length = len;
      if (type == ATTRINT) {
        if (len == 1) {
          attr[n].vals.integer  = ival;
        } else {
          attr[n].vals.integers = ivec;
        }
      } else if ((type == ATTRREAL) || (type == ATTRCSYS)) {
        if (len == 1) {
          attr[n].vals.real  = rval;
        } else {
          attr[n].vals.reals = rvec;
        }
      } else {
        attr[n].vals.string = string;
      }
      n++;
    }
  }

  if (attrs != NULL) {
    attrs->nattrs = n;
    attrs->attrs  = attr;
    attrs->nseqs  = 0;
    attrs->seqs   = NULL;
    if (nseq != 0) EG_attrBuildSeq(attrs);
    obj->attrs    = attrs;
  }
}


void
EG_initOCC()
{
  int  np;
  char *env;

  env = getenv("EMPnumProc");
  np  = 2;
  if (env != NULL) np = atoi(env);
  if (np > 1) BOPAlgo_Options::SetParallelMode(Standard_True);

#ifdef OCC_VERSION_DEVELOPMENT
  if (strcmp("ESP", OCC_VERSION_DEVELOPMENT) == 0) return;
#endif
  if (((OCC_VERSION_MAJOR       != 7) || (OCC_VERSION_MINOR != 3) ||
       (OCC_VERSION_MAINTENANCE == 0)) &&
      ((OCC_VERSION_MAJOR       != 7) || (OCC_VERSION_MINOR != 4) ||
       (OCC_VERSION_MAINTENANCE == 0)) &&
      ((OCC_VERSION_MAJOR       != 7) || (OCC_VERSION_MINOR != 6) ||
       (OCC_VERSION_MAINTENANCE == 0)))
    printf(" EGADS WARNING: OpenCASCADE %d.%d.%d NOT an Authorized Release!\n",
           OCC_VERSION_MAJOR, OCC_VERSION_MINOR, OCC_VERSION_MAINTENANCE);
}


static int
EG_readTess(FILE *fp, egObject *body, egObject **tess)
{
  int     i, j, ir, status, nnode, nedge, nface, n[3], len, ntri, nattr;
  int     *ptype, *pindex, *tris, *tric, nFace;
  double  *xyz, *param;
  egObject *obj;

  *tess  = NULL;
  status = EG_getBodyTopos(body, NULL, NODE, &nnode, NULL);
  if (status != EGADS_SUCCESS) return status;
  if (body->oclass == EBODY) {
    status = EG_getBodyTopos(body, NULL, EEDGE, &nedge, NULL);
    if (status != EGADS_SUCCESS) return status;
    status = EG_getBodyTopos(body, NULL, EFACE, &nface, NULL);
    if (status != EGADS_SUCCESS) return status;
  } else {
    status = EG_getBodyTopos(body, NULL, EDGE, &nedge, NULL);
    if (status != EGADS_SUCCESS) return status;
    status = EG_getBodyTopos(body, NULL, FACE, &nface, NULL);
    if (status != EGADS_SUCCESS) return status;
  }

  ir = fscanf(fp, "%d %d %d", &n[0], &n[1], &n[2]);
  if (ir != 3) {
    printf(" EGADS Error: Header with only %d words (EG_readTess)!\n", ir);
    return EGADS_INDEXERR;
  }
  if ((nnode != n[0]) || (nedge != n[1]) || (nface != n[2])) {
    printf(" EGADS Error: Count mismatch %d %d  %d %d  %d %d (EG_readTess)!\n",
           nnode, n[0], nedge, n[1], nface, n[2]);
    return EGADS_INDEXERR;
  }

  /* initialize the Tessellation Object */
  status = EG_initTessBody(body, tess);
  if (status != EGADS_SUCCESS) return status;
  EG_dereferenceTopObj(body, *tess);

  /* do the Edges */
  for (i = 0; i < nedge; i++) {
    len = 0;
    fscanf(fp, "%d", &len);
    if (len == 0) continue;
    xyz   = (double *) malloc(3*len*sizeof(double));
    param = (double *) malloc(  len*sizeof(double));
    if ((xyz == NULL) || (param == NULL)) {
      printf(" EGADS Error: malloc on Edge %d -- len = %d (EG_readTess)!\n",
               i+1, len);
      if (xyz   != NULL) free(xyz);
      if (param != NULL) free(param);
      EG_deleteObject(*tess);
      *tess = NULL;
      return EGADS_MALLOC;
    }
    for (j = 0; j < len; j++) {
      ir = fscanf(fp, "%le %le %le %le", &xyz[3*j  ], &xyz[3*j+1], &xyz[3*j+2],
                  &param[j]);
      if (ir != 4) {
        printf(" EGADS Error: %d/%d Read got %d out of 4 (EG_readTess)!\n",
               j+1, len, ir);
        free(xyz);
        free(param);
        EG_deleteObject(*tess);
        *tess = NULL;
        return EGADS_READERR;
      }
    }
    status = EG_setTessEdge(*tess, i+1, len, xyz, param);
    free(xyz);
    free(param);
    if (status != EGADS_SUCCESS) {
      printf(" EGADS Error: EG_setTessEdge %d = %d (EG_readTess)!\n",
             i+1, status);
      EG_deleteObject(*tess);
      *tess = NULL;
      return status;
    }
  }

  /* do the Faces */
  for (nFace = i = 0; i < nface; i++) {
    len = ntri = 0;
    fscanf(fp, "%d %d", &len, &ntri);
    if ((len == 0) || (ntri == 0)) {
      continue;
    }
    nFace++; /* Count total number of faces with tessellations */
    xyz    = (double *) malloc(3*len*sizeof(double));
    param  = (double *) malloc(2*len*sizeof(double));
    ptype  = (int *)    malloc(  len* sizeof(int));
    pindex = (int *)    malloc(  len* sizeof(int));
    tris   = (int *)    malloc(3*ntri*sizeof(int));
    tric   = (int *)    malloc(3*ntri*sizeof(int));
    if ((xyz    == NULL) || (param == NULL) || (ptype == NULL) ||
        (pindex == NULL) || (tris  == NULL) || (tric  == NULL)) {
      printf(" EGADS Error: malloc on Face %d -- lens = %d %d (EG_readTess)!\n",
             i+1, len, ntri);
      if (xyz    != NULL) free(xyz);
      if (param  != NULL) free(param);
      if (ptype  != NULL) free(ptype);
      if (pindex != NULL) free(pindex);
      if (tris   != NULL) free(tris);
      if (tric   != NULL) free(tric);
      EG_deleteObject(*tess);
      *tess = NULL;
      return EGADS_MALLOC;
    }
    for (j = 0; j < len; j++) {
      ir = fscanf(fp, "%le %le %le %le %le %d %d", &xyz[3*j], &xyz[3*j+1],
                  &xyz[3*j+2], &param[2*j], &param[2*j+1], &ptype[j],
                  &pindex[j]);
      if (ir != 7) {
        printf(" EGADS Error: %d/%d Read got %d out of 7 (EG_readTess)!\n",
               j+1, len, ir);
        free(xyz);
        free(param);
        free(ptype);
        free(pindex);
        free(tris);
        free(tric);
        EG_deleteObject(*tess);
        *tess = NULL;
        return EGADS_READERR;
      }
    }
    for (j = 0; j < ntri; j++) {
      ir = fscanf(fp, "%d %d %d %d %d %d", &tris[3*j], &tris[3*j+1],
                  &tris[3*j+2], &tric[3*j], &tric[3*j+1], &tric[3*j+2]);
      if (ir != 6) {
        printf(" EGADS Error: %d/%d Read got %d out of 6 (EG_readTess)!\n",
               j+1, len, ir);
        free(xyz);
        free(param);
        free(ptype);
        free(pindex);
        free(tris);
        free(tric);
        EG_deleteObject(*tess);
        *tess = NULL;
        return EGADS_READERR;
      }
    }
    status = EG_setTessFace(*tess, i+1, len, xyz, param, ntri, tris);
    free(xyz);
    free(param);
    free(ptype);
    free(pindex);
    free(tris);
    free(tric);
    if (status != EGADS_SUCCESS) {
      printf(" EGADS Warning: EG_setTessFace %d = %d (EG_readTess)!\n",
             i+1, status);
/*    EG_deleteObject(*tess);
      *tess = NULL;
      return status; */
    }
  }

  if (nFace == 0) {
    /* Special case where there are no face tessellations */
    egTessel *btess = (egTessel *) (*tess)->blind;
    btess->nFace = 0;
    EG_free(btess->tess2d);
    btess->tess2d = NULL;
  }

  /* close up the open tessellation */
  status = EG_statusTessBody(*tess, &obj, &i, &len);
  if (status == EGADS_OUTSIDE) {
    printf(" EGADS Warning: Tessellation Object is incomplete (EG_readTess)!\n");
    egTessel *btess = (egTessel *) (*tess)->blind;
    for (i = 0; i < btess->nFace; i++) {
      if ((btess->tess2d[i].npts == 0) || (btess->tess2d[i].ntris == 0)) {
        printf("                Face %d -- lens = %d %d (EG_readTess)!\n",
               i+1, btess->tess2d[i].npts, btess->tess2d[i].ntris);
      }
    }
    btess->done = 1;
  } else if (status != EGADS_SUCCESS) {
    printf(" EGADS Error: EG_statusTessBody = %d (EG_readTess)!\n", status);
    EG_deleteObject(*tess);
    *tess = NULL;
    return status;
  }
  if ((status != EGADS_OUTSIDE) && (i != 1)) {
    printf(" EGADS Warning: Tessellation Object is %d (EG_readTess)!\n", i);
/*  EG_deleteObject(*tess);
    *tess = NULL;
    return EGADS_TESSTATE;  */
    egTessel *btess = (egTessel *) (*tess)->blind;
    btess->done = 1;
  }

  /* attach the attributes */
  fscanf(fp, "%d\n", &nattr);
  if (nattr != 0) EG_readAttrs(*tess, nattr, fp);

  return EGADS_SUCCESS;
}


struct MatchShape
{
  typedef NCollection_IndexedDataMap<Handle(Geom_Curve),    TopoDS_Shape>
                                    GeomCurve_IndexedDataMap;
  typedef NCollection_IndexedDataMap<Handle(Geom_Surface),  TopoDS_Shape>
                                    GeomSurface_IndexedDataMap;
  typedef NCollection_IndexedDataMap<Handle(TopoDS_TShape), TopoDS_Shape>
                                    TShape_IndexedDataMap;

  MatchShape(const TopoDS_Shape& aShape)
  {
    TopExp::MapShapes(aShape, shapes);

    for (int i = 1; i <= shapes.Extent(); i++) {
      const TopoDS_Shape& shape = shapes(i);
      tshapes.Add(shape.TShape(), shape);
      if (shape.ShapeType() == TopAbs_EDGE) {
        TopLoc_Location    loc;
        Standard_Real      first, last;
        Handle(Geom_Curve) Crv = BRep_Tool::Curve(TopoDS::Edge(shape), loc,
                                                  first, last);
        if (!Crv.IsNull())
          curves.Add(Crv, shape);
      } else if (shape.ShapeType() == TopAbs_FACE) {
        Handle(Geom_Surface) Srf = BRep_Tool::Surface(TopoDS::Face(shape));
        if (!Srf.IsNull())
          surfaces.Add(Srf, shape);
      }
    }
  }

  void Find(TopoDS_Shape& S) const
  {
    if (shapes.Contains(S)) return;

    Standard_Boolean found = false;
    if (S.ShapeType() == TopAbs_EDGE) {

      // Find Edge by matching underlying curve
      TopoDS_Edge edge = TopoDS::Edge (S);

      TopLoc_Location    loc;
      Standard_Real      first, last;
      Handle(Geom_Curve) Crv = BRep_Tool::Curve(edge, loc, first, last);
      if (!Crv.IsNull()) {
        found = curves.FindFromKey(Crv, S);
      }

    } else if (S.ShapeType() == TopAbs_WIRE) {

      // Find Edge by matching underlying curve (yes it's in a wire...)
      TopTools_IndexedMapOfShape edges;
      TopExp::MapShapes(S, TopAbs_EDGE, edges);

      if (edges.Extent() == 1) {
        TopoDS_Edge edge = TopoDS::Edge (edges(1));

        TopLoc_Location    loc;
        Standard_Real      first, last;
        Handle(Geom_Curve) Crv = BRep_Tool::Curve(edge, loc, first, last);
        if (!Crv.IsNull()) {
          found = curves.FindFromKey(Crv, S);
        }
      }

    } else if (S.ShapeType() == TopAbs_FACE) {

      // Find Face by matching underlying surface
      Handle(Geom_Surface) Srf = BRep_Tool::Surface(TopoDS::Face(S));

      if (!Srf.IsNull()) {
        found = surfaces.FindFromKey(Srf, S);
      }
    }

    // find the IsPartner shape if necessary
    if (!found) {
      tshapes.FindFromKey(S.TShape(), S);
    }
  }

protected:

  TopTools_IndexedMapOfShape shapes;
  TShape_IndexedDataMap tshapes;
  GeomCurve_IndexedDataMap curves;
  GeomSurface_IndexedDataMap surfaces;
};


// Taken from
// XSControl_TransferReader::EntityFromShapeResult
// XSControl_TransferReader::EntitiesFromShapeList
static void
EG_labelsFromTransferReader
       (const TopoDS_Shape& aShape,
        const Handle(XSControl_TransferReader)& TR,
        egadsShapeData& labels)
{
  const Handle(Transfer_TransientProcess)& TP = TR->TransientProcess();
  const Handle(Interface_InterfaceModel)& Model = TR->Model();

  MatchShape match(aShape);

  Standard_Integer i, j, nb;

  XSControl_Utils xu;
  if (!TP.IsNull()) {
    nb = TP->NbMapped();
    for (j = 1; j <= nb; j++) {
      Handle(Standard_Transient) ent = TP->Mapped(j);

      Handle(TCollection_HAsciiString) name;

      // First try STEP file
      const Handle(StepRepr_RepresentationItem) aReprItem =
            Handle(StepRepr_RepresentationItem)::DownCast(ent);
      if (!aReprItem.IsNull()) {
        name = aReprItem->Name();
        if (name.IsNull() || name->IsEmpty()) {
          // Check underlying geometry
          Handle(StepShape_VertexPoint)
            vertPoint = Handle(StepShape_VertexPoint)::DownCast(aReprItem);
          if (!vertPoint.IsNull())
            name = vertPoint->VertexGeometry()->Name();
          else {
            Handle(StepShape_EdgeCurve)
              edgeCurve = Handle(StepShape_EdgeCurve)::DownCast(aReprItem);
            if (!edgeCurve.IsNull())
              name = edgeCurve->EdgeGeometry()->Name();
            else {
              Handle(StepShape_FaceSurface)
                faceSurf = Handle(StepShape_FaceSurface)::DownCast(aReprItem);
              if (!faceSurf.IsNull())
                name = faceSurf->FaceGeometry()->Name();
            }
          }
        }
      } else {
        // Try IGES file

#if CASVER < 770
        // check each curve for an IGESSolid_EdgeList
        const Handle(IGESSolid_EdgeList)& edgeList =
              Handle(IGESSolid_EdgeList)::DownCast(ent);
        if (!edgeList.IsNull()) {
          const Handle(TransferBRep_ShapeListBinder)& binder =
                Handle(TransferBRep_ShapeListBinder)::DownCast(TP->Find(ent));
          if (binder.IsNull()) continue;
          if (binder->NbShapes() != edgeList->NbEdges()) continue;

          for (Standard_Integer iedge = 1; iedge <= edgeList->NbEdges(); iedge++) {
            Handle(IGESData_IGESEntity) Curve = edgeList->Curve(iedge);

            name.Nullify();
            if (!Curve.IsNull() && Curve->HasName()) name = Curve->NameValue();
            if (name.IsNull() || name->IsEmpty()) continue;

            TopoDS_Shape sh = binder->Edge(iedge);
            if (sh.IsNull()) continue;
            match.Find(sh);
            labels.Add(sh, name->ToCString());
          }
          name.Nullify();
        }
        
        // check each curve for an IGESGeom_CompositeCurve
        const Handle(IGESGeom_CompositeCurve)& compCurve =
              Handle(IGESGeom_CompositeCurve)::DownCast(ent);
        if (!compCurve.IsNull()) {
          const Handle(TransferBRep_ShapeBinder)& binder =
                Handle(TransferBRep_ShapeBinder)::DownCast(TP->Find(ent));
          if (binder.IsNull()) continue;

          TopoDS_Shape wire = binder->Result();
          if (wire.IsNull()) continue;
          match.Find(wire);

          TopExp_Explorer Exp;
          Exp.Init(wire, TopAbs_EDGE);

          for (Standard_Integer iedge = 1;
              iedge <= compCurve->NbCurves(); iedge++, Exp.Next()) {
            Handle(IGESData_IGESEntity) Curve = compCurve->Curve(iedge);

            name.Nullify();
            if (!Curve.IsNull() && Curve->HasName()) name = Curve->NameValue();
            if (name.IsNull() || name->IsEmpty()) continue;

            TopoDS_Shape sh = Exp.Current();
            match.Find(sh);
            labels.Add(sh, name->ToCString());
          }
          name.Nullify();
        }
#endif

        const Handle(IGESData_IGESEntity)& shapeEntity =
              Handle(IGESData_IGESEntity)::DownCast(ent);
        if (!shapeEntity.IsNull() && shapeEntity->HasName())
          name = shapeEntity->NameValue();
      }
      if (name.IsNull() || name->IsEmpty()) continue;

      TopoDS_Shape sh = TransferBRep::ShapeResult (TP,ent);
      if (sh.IsNull()) continue;

      match.Find(sh);

      labels.Add(sh, name->ToCString());
      /* printf(" Name found: %s\n",
              TopAbs::ShapeTypeToString(sh.ShapeType()), name ); */
    }
  }

  if (Model.IsNull()) return;

  nb = Model->NbEntities();
  for (i = 1; i <= nb; i ++) {
    Handle(Transfer_ResultFromModel) rec = TR->ResultFromNumber (i);
    if (rec.IsNull()) continue;

    Handle(TColStd_HSequenceOfTransient) list = rec->Results (2);
    Standard_Integer ir,nr = list->Length();
    for (ir = 1; ir <= nr; ir ++) {
      DeclareAndCast(Transfer_ResultFromTransient,rft,list->Value(ir));
      if (rft.IsNull()) continue;
      Handle(Standard_Transient) ent = rft->Start();
      if (ent.IsNull()) continue;

      Handle(TCollection_HAsciiString) name;

      // First try STEP file
      Handle(StepRepr_RepresentationItem) aReprItem;
      aReprItem = Handle(StepRepr_RepresentationItem)::DownCast(ent);
      if (!aReprItem.IsNull()) {
        name = aReprItem->Name();
        if (name.IsNull() || name->IsEmpty()) {
          // Check underlying geometry
          Handle(StepShape_VertexPoint)
            vertPoint = Handle(StepShape_VertexPoint)::DownCast(aReprItem);
          if (!vertPoint.IsNull())
            name = vertPoint->VertexGeometry()->Name();
          else {
            Handle(StepShape_EdgeCurve)
              edgeCurve = Handle(StepShape_EdgeCurve)::DownCast(aReprItem);
            if (!edgeCurve.IsNull())
              name = edgeCurve->EdgeGeometry()->Name();
            else {
              Handle(StepShape_FaceSurface)
                faceSurf = Handle(StepShape_FaceSurface)::DownCast(aReprItem);
              if (!faceSurf.IsNull())
                name = faceSurf->FaceGeometry()->Name();
            }
          }
        }
      } else {
        // Try IGES file

#if CASVER < 770
        // check each curve for an IGESSolid_EdgeList
        const Handle(IGESSolid_EdgeList)& edgeList =
              Handle(IGESSolid_EdgeList)::DownCast(ent);
        if (!edgeList.IsNull()) {
          const Handle(TransferBRep_ShapeListBinder)& binder =
                Handle(TransferBRep_ShapeListBinder)::DownCast(TP->Find(ent));
          if (binder.IsNull()) continue;
          if (binder->NbShapes() != edgeList->NbEdges()) continue;

          for (Standard_Integer iedge = 1; iedge <= edgeList->NbEdges(); iedge++) {
            Handle(IGESData_IGESEntity) Curve = edgeList->Curve(iedge);

            name.Nullify();
            if (!Curve.IsNull() && Curve->HasName()) name = Curve->NameValue();
            if (name.IsNull() || name->IsEmpty()) continue;

            TopoDS_Shape sh = binder->Edge(iedge);
            if (sh.IsNull()) continue;
            match.Find(sh);
            labels.Add(sh, name->ToCString());
          }
          name.Nullify();
        }

        // check each curve for an IGESGeom_CompositeCurve
        const Handle(IGESGeom_CompositeCurve)& compCurve =
              Handle(IGESGeom_CompositeCurve)::DownCast(ent);
        if (!compCurve.IsNull()) {
          const Handle(TransferBRep_ShapeBinder)& binder =
                Handle(TransferBRep_ShapeBinder)::DownCast(TP->Find(ent));
          if (binder.IsNull()) continue;

          TopoDS_Shape wire = binder->Result();
          if (wire.IsNull()) continue;
          match.Find(wire);

          TopExp_Explorer Exp;
          Exp.Init(wire, TopAbs_EDGE);

          for (Standard_Integer iedge = 1;
              iedge <= compCurve->NbCurves(); iedge++, Exp.Next()) {
            Handle(IGESData_IGESEntity) Curve = compCurve->Curve(iedge);

            name.Nullify();
            if (!Curve.IsNull() && Curve->HasName()) name = Curve->NameValue();
            if (name.IsNull() || name->IsEmpty()) continue;

            TopoDS_Shape sh = Exp.Current();
            match.Find(sh);
            labels.Add(sh, name->ToCString());
          }
          name.Nullify();
        }
#endif

        const Handle(IGESData_IGESEntity)& shapeEntity =
              Handle(IGESData_IGESEntity)::DownCast(ent);
        if (!shapeEntity.IsNull() && shapeEntity->HasName())
          name = shapeEntity->NameValue();
      }
      if (name.IsNull() || name->IsEmpty()) continue;

      TopoDS_Shape sh = xu.BinderShape (rft->Binder());
      if (sh.IsNull()) continue;
      
      match.Find(sh);

      labels.Add(sh, name->ToCString());
      /* printf(" Name found: %s\n",
              TopAbs::ShapeTypeToString(sh.ShapeType()), name ); */
    }
  }
}

// Based on STEPConstruct_Styles::DecodeColor
Standard_Boolean
EG_DecodeColor (const Handle(StepVisual_Colour) &Colour, Quantity_Color &Col)
{
  if ( Colour->IsKind(STANDARD_TYPE(StepVisual_ColourRgb)) ) {
    Handle(StepVisual_ColourRgb)
      rgb = Handle(StepVisual_ColourRgb)::DownCast(Colour);
    if (rgb->Red() > 1.0 || rgb->Green() > 1.0 || rgb->Blue() > 1.0) {
      Standard_Real norm = rgb->Red();
      if(norm<rgb->Green()) norm = rgb->Green();
      if(norm<rgb->Blue()) norm = rgb->Blue();
      Col.SetValues(rgb->Red()/norm, rgb->Green()/norm,
                    rgb->Blue()/norm, Quantity_TOC_RGB);
    }
    else
      Col.SetValues(rgb->Red(), rgb->Green(), rgb->Blue(), Quantity_TOC_RGB);
    return Standard_True;
  }
  else if ( Colour->IsKind (STANDARD_TYPE(StepVisual_PreDefinedColour)) ) {
    Handle(StepVisual_PreDefinedColour) pdc =
      Handle(StepVisual_PreDefinedColour)::DownCast ( Colour );
    Handle(StepVisual_PreDefinedItem) pdi = pdc->GetPreDefinedItem();
    const TCollection_AsciiString name = pdi->Name()->String();
    if      ( name.IsEqual ( "red"     ) ) Col.SetValues ( Quantity_NOC_RED );
    else if ( name.IsEqual ( "green"   ) ) Col.SetValues ( Quantity_NOC_GREEN );
    else if ( name.IsEqual ( "blue"    ) ) Col.SetValues ( Quantity_NOC_BLUE1 );
    else if ( name.IsEqual ( "yellow"  ) ) Col.SetValues ( Quantity_NOC_YELLOW );
    else if ( name.IsEqual ( "magenta" ) ) Col.SetValues ( Quantity_NOC_MAGENTA1 );
    else if ( name.IsEqual ( "cyan"    ) ) Col.SetValues ( Quantity_NOC_CYAN1 );
    else if ( name.IsEqual ( "black"   ) ) Col.SetValues ( Quantity_NOC_BLACK );
    else if ( name.IsEqual ( "white"   ) ) Col.SetValues ( Quantity_NOC_WHITE );
    else if ( name.IsEqual ( "lred"    ) ) Col.SetValues(1.0, 0.5, 0.5, Quantity_TOC_RGB);
    else if ( name.IsEqual ( "lgreen"  ) ) Col.SetValues(0.5, 1.0, 0.5, Quantity_TOC_RGB);
    else if ( name.IsEqual ( "lblue"   ) ) Col.SetValues(0.5, 0.5, 1.0, Quantity_TOC_RGB);
    else {
#ifdef OCCT_DEBUG
      std::cout << "Error: color name \"" << name << "\" is not recognized" <<
      std::endl;
#endif
      return Standard_False;
    }
    return Standard_True;
  }
  return Standard_False;
}


#ifdef SETP_COLOR_NOTNEEDED
// Based on STEPCAFControl_Reader SetAssemblyComponentStyle
static void
EG_setAssemblyComponentStyle(const Handle(Transfer_TransientProcess) &theTP,
                             const STEPConstruct_Styles& theStyles,
                             const Handle(StepVisual_ContextDependentOverRidingStyledItem)& theStyle,
                             egadsShapeData& shapeData)
{
  if (theStyle.IsNull()) return;

  Handle(StepVisual_Colour) aSurfCol, aBoundCol, aCurveCol, aRenderCol;
  Standard_Real aRenderTransp;
  // check if it is component style
  Standard_Boolean anIsComponent = Standard_False;
  if (!theStyles.GetColors(theStyle, aSurfCol, aBoundCol, aCurveCol, aRenderCol,
                           aRenderTransp, anIsComponent))
    return;

  TopLoc_Location aLoc; // init;
  // find shape
  TopoDS_Shape aShape;
  Handle(Transfer_Binder) aBinder;
  Handle(StepShape_ShapeRepresentation)
    aRepr = Handle(StepShape_ShapeRepresentation)::DownCast(theStyle->ItemAP242().Value());
  if (aRepr.IsNull())
    return;
  Handle(StepRepr_ShapeRepresentationRelationship) aSRR;
  Interface_EntityIterator aSubs = theTP->Graph().Sharings(aRepr);
  for (aSubs.Start(); aSubs.More(); aSubs.Next()) {
    const Handle(Standard_Transient)& aSubsVal = aSubs.Value();
    if (aSubsVal->IsKind (STANDARD_TYPE(StepRepr_ShapeRepresentationRelationship)))
    {
      // NB: C cast is used instead of DownCast() to improve performance on some cases.
      // This saves ~10% of elapsed time on "testgrid perf de bug29* -parallel 0".
      aSRR = (StepRepr_ShapeRepresentationRelationship*)(aSubsVal.get());
    }
  }

  aBinder = theTP->Find(aSRR);
  if (!aBinder.IsNull() ) {
    aShape = TransferBRep::ShapeResult (aBinder);
  }
  if (aShape.IsNull())
    return;

  if(!aSurfCol.IsNull() || !aBoundCol.IsNull() || !aCurveCol.IsNull() ||
     !aRenderCol.IsNull()) {
    Quantity_Color aSCol,aBCol,aCCol,aRCol;
    if(!aSurfCol.IsNull()) {
      theStyles.DecodeColor(aSurfCol,aSCol);
      if (aShape.ShapeType() == TopAbs_FACE) {
        shapeData.Add(aShape, aSCol);
      }
    }
    if(!aBoundCol.IsNull())
      theStyles.DecodeColor(aBoundCol,aBCol);
    if(!aCurveCol.IsNull()) {
      theStyles.DecodeColor(aCurveCol,aCCol);
      if (aShape.ShapeType() == TopAbs_EDGE) {
        shapeData.Add(aShape, aCCol);
      }
    }
    if(!aRenderCol.IsNull()) {
      theStyles.DecodeColor(aRenderCol,aRCol);
    }
  }
}
#endif


// Based on STEPCAFControl_Reader SetStyle
static void
EG_setStyle(const MatchShape& match, const Handle(XSControl_WorkSession) &theWS,
            const STEPConstruct_Styles& theStyles,
            const Handle(StepVisual_StyledItem)& theStyle,
            egadsShapeData& shapeData)
{
  const Handle(Transfer_TransientProcess)
        &aTP = theWS->TransferReader()->TransientProcess();
#ifdef SETP_COLOR_NOTNEEDED
  if (Handle(StepVisual_OverRidingStyledItem)
        anOverridingStyle = Handle(StepVisual_OverRidingStyledItem)::DownCast(theStyle)) {
    EG_setStyle(theWS, theStyles, anOverridingStyle->OverRiddenStyle(), shapeData);
    if (Handle(StepVisual_ContextDependentOverRidingStyledItem) anAssemblyComponentStyle =
        Handle(StepVisual_ContextDependentOverRidingStyledItem)::DownCast(theStyle)) {
      EG_setAssemblyComponentStyle(aTP, theStyles, anAssemblyComponentStyle, shapeData);
      return;
    }
  }
#endif

  Handle(StepVisual_Colour) aSurfCol, aBoundCol, aCurveCol;
  // check if it is component style
  Standard_Boolean anIsComponent = Standard_False;
#if CASVER >= 760
  Handle(StepVisual_Colour) aRenderCol;
  Standard_Real aRenderTransp;
  if (!theStyles.GetColors(theStyle, aSurfCol, aBoundCol, aCurveCol, aRenderCol,
                           aRenderTransp, anIsComponent)) return;
#else
  if (!theStyles.GetColors(theStyle, aSurfCol, aBoundCol, aCurveCol,
                           anIsComponent)) return;
#endif

  // collect styled items
  NCollection_Vector<StepVisual_StyledItemTarget> anItems;
  if (!theStyle->ItemAP242().IsNull()) {
    anItems.Append(theStyle->ItemAP242());
  }

  for (Standard_Integer itemIt = 0; itemIt < anItems.Length(); itemIt++) {
    Standard_Integer anIndex = aTP->MapIndex(anItems.Value(itemIt).Value());
    TopoDS_Shape aShape;
    if (anIndex > 0) {
      Handle(Transfer_Binder) aBinder = aTP->MapItem(anIndex);
      aShape = TransferBRep::ShapeResult(aBinder);
    }
    if (aShape.IsNull()) continue;

    match.Find(aShape);

    Quantity_Color aCCol;
    if (!aCurveCol.IsNull()) {
      EG_DecodeColor(aCurveCol, aCCol);

      if (aShape.ShapeType() == TopAbs_EDGE) {

        shapeData.Add(aShape, aCCol);

      } else if (aShape.ShapeType() == TopAbs_FACE) {
        // Apply Face Curve colors to Edges
        TopExp_Explorer Exp;
        for (Exp.Init(aShape, TopAbs_EDGE); Exp.More(); Exp.Next()) {
          TopoDS_Shape sh = Exp.Current();
          match.Find(sh);
          shapeData.Add(sh, aCCol);
        }
      }
    }

    Quantity_Color aSCol;
    if (!aSurfCol.IsNull()) {
      EG_DecodeColor(aSurfCol, aSCol);
    }

    if (!aSurfCol.IsNull() && aShape.ShapeType() == TopAbs_FACE) {
      shapeData.Add(aShape, aSCol);
    }

    if ((!aCurveCol.IsNull() || !aSurfCol.IsNull()) &&
        (aShape.ShapeType() == TopAbs_SHELL ||
         aShape.ShapeType() == TopAbs_SOLID ||
         aShape.ShapeType() == TopAbs_COMPOUND)) { // WIRE bodies come in as COMPOUND

      Quantity_Color  *CCol = NULL, *SCol = NULL;
      if (!aCurveCol.IsNull()) CCol = new Quantity_Color(aCCol);
      if (!aSurfCol.IsNull() ) SCol = new Quantity_Color(aSCol);

      shapeData.Add(aShape, CCol, SCol);
    }

//    Quantity_Color aBCol;
//    if (!aBoundCol.IsNull())
//      theStyles.DecodeColor(aBoundCol, aBCol);
//    if (!aRenderCol.IsNull())
//      theStyles.DecodeColor(aRenderCol, aRCol);
  }
}


// Based on STEPCAFControl_Reader IsOverriden
static Standard_Boolean
EG_stepIsOverriden(const Interface_Graph& theGraph,
                   const Handle(StepVisual_StyledItem)& theStyle,
                   Standard_Boolean theIsRoot)
{
  Interface_EntityIterator aSubs = theGraph.Sharings (theStyle);
  aSubs.Start();
  for(; aSubs.More(); aSubs.Next()) {
    Handle(StepVisual_OverRidingStyledItem)
      anOverRidingStyle = Handle(StepVisual_OverRidingStyledItem)::DownCast(aSubs.Value());
    if(!anOverRidingStyle.IsNull()) {
      if(!theIsRoot) {
        return Standard_True;
      }
      // for root style returns true only if it is overridden by other root style
      const Handle(Standard_Transient) anItem = anOverRidingStyle->ItemAP242().Value();
      if(!anItem.IsNull() &&
         anItem->IsKind(STANDARD_TYPE(StepShape_ShapeRepresentation))) {
        return Standard_True;
      }
    }
  }
  return Standard_False;
}


// Based on STEPCAFControl_Reader::ReadColors
static void
EG_stepColors(const TopoDS_Shape& aShape,
              const Handle(XSControl_WorkSession)& WS,
              egadsShapeData& shapeData)
{

  STEPConstruct_Styles Styles(WS);
  if (!Styles.LoadStyles()) return;

  MatchShape match(aShape);

  const Interface_Graph& aGraph = Styles.Graph ();

  Standard_Boolean anIsRootStyle = Standard_True;

  Standard_Integer nb;
#if CASVER >= 780
  // parse and search for color attributes
  nb = Styles.NbRootStyles();
  // apply root styles earlier, as they can be overridden
  // function IsOverriden for root style returns true only if it is
  //          overridden by other root style
  for(Standard_Integer i = 1; i <= nb; i++)
  {
    Handle(StepVisual_StyledItem) Style = Styles.RootStyle(i);
    // check that style is overridden by other root style
    if (!EG_stepIsOverriden (aGraph, Style, anIsRootStyle)) {
      EG_setStyle (match, WS, Styles, Style, shapeData);
    }
  }
#endif

  nb = Styles.NbStyles();
  // apply leaf styles, they can override root styles
  anIsRootStyle = Standard_False;
  for(Standard_Integer i = 1; i <= nb; i++)
  {
    Handle(StepVisual_StyledItem) Style = Styles.Style(i);
    // check that style is overridden
    if (!EG_stepIsOverriden(aGraph, Style, anIsRootStyle)) {
      EG_setStyle(match, WS, Styles, Style, shapeData);
    }
  }

}


// Based on IGESCAFControl::DecodeColor
Quantity_Color
EG_igesDecodeColor(const Standard_Integer color)
{
  switch (color) {
    case 1:  return Quantity_Color ( Quantity_NOC_BLACK );
    case 2:  return Quantity_Color ( Quantity_NOC_RED );
    case 3:  return Quantity_Color ( Quantity_NOC_GREEN );
    case 4:  return Quantity_Color ( Quantity_NOC_BLUE1 );
    case 5:  return Quantity_Color ( Quantity_NOC_YELLOW );
    case 6:  return Quantity_Color ( Quantity_NOC_MAGENTA1 );
    case 7:  return Quantity_Color ( Quantity_NOC_CYAN1 );
    case 8:
    default: return Quantity_Color ( Quantity_NOC_WHITE );
  }
}


// Based on IGESCAFControl::EncodeColor
Standard_Integer
EG_igesEncodeColor(const Quantity_Color &col)
{
  Standard_Integer code = 0;
  if      (Abs(col.Red()   - 1.0) <= col.Epsilon()) code |= 0x001;
  else if (Abs(col.Red())         >  col.Epsilon()) return 0;
  if      (Abs(col.Green() - 1.0) <= col.Epsilon()) code |= 0x010;
  else if (Abs(col.Green())       >  col.Epsilon()) return 0;
  if      (Abs(col.Blue()  - 1.0) <= col.Epsilon()) code |= 0x100;
  else if (Abs(col.Blue())        >  col.Epsilon()) return 0;

  switch (code) {
    case 0x000: return 1;
    case 0x001: return 2;
    case 0x010: return 3;
    case 0x100: return 4;
    case 0x011: return 5;
    case 0x101: return 6;
    case 0x110: return 7;
    case 0x111:
    default   : return 8;
  }
}


static void
EG_igesCheckColorRange(Standard_Real& theCol)
{
  if (theCol <   0.0) theCol =   0.0;
  if (theCol > 100.0) theCol = 100.0;
}


// EG_igesEntityColor & EG_igesColors based on IGESCAFControl_Reader::Transfer
static void
EG_igesEntityColor(const Handle(IGESData_IGESEntity)& ent,
                   const MatchShape& match,
                   const Handle(Transfer_TransientProcess) &TP,
                   TopoDS_Shape S, egadsShapeData& shapeData)
{
  Standard_Boolean IsColor = Standard_False;
  Quantity_Color col;
  // read colors
  if(ent->DefColor() == IGESData_DefValue ||
     ent->DefColor() == IGESData_DefReference) {
    // color is assigned
    // decode color and set to document
    IsColor = Standard_True;
    if (ent->DefColor() == IGESData_DefValue) {
      col = EG_igesDecodeColor(ent->RankColor());
    } else {
      Handle(IGESGraph_Color)
        color = Handle(IGESGraph_Color)::DownCast(ent->Color());
      if (color.IsNull()) {
        IsColor = Standard_False;
      } else {
        Standard_Real r, g, b;
        color->RGBIntensity(r, g, b);
        EG_igesCheckColorRange(r);
        EG_igesCheckColorRange(g);
        EG_igesCheckColorRange(b);
        col.SetValues(0.01*r, 0.01*g, 0.01*b, Quantity_TOC_RGB);
      }
    }
  }
  if (!IsColor) return;

  // Find the shape
  if (S.IsNull()) {
    Handle(Transfer_Binder) binder = TP->Find(ent);
    if (binder.IsNull()) return;
    S = TransferBRep::ShapeResult(binder);
    if (S.IsNull()) return;
  }

  match.Find(S);

  shapeData.Add(S, col);
}


static void
EG_igesColors(const TopoDS_Shape& aShape,
              const Handle(XSControl_WorkSession)& WS,
              egadsShapeData& shapeData)
{
  MatchShape match(aShape);

  Handle(IGESData_IGESModel)
    aModel = Handle(IGESData_IGESModel)::DownCast(WS->Model());
  const Handle(XSControl_TransferReader)  &TR = WS->TransferReader();
  const Handle(Transfer_TransientProcess) &TP = TR->TransientProcess();

  Standard_Integer nb = aModel->NbEntities();
  for(Standard_Integer i = 1; i <= nb; i++) {
    Handle(IGESData_IGESEntity)
      ent = Handle(IGESData_IGESEntity)::DownCast(aModel->Value(i));
    if ( ent.IsNull() ) continue;

    TopoDS_Shape S;
    EG_igesEntityColor(ent, match, TP, S, shapeData);

    // check each curve for an IGESSolid_EdgeList
    const Handle(IGESSolid_EdgeList)& edgeList =
          Handle(IGESSolid_EdgeList)::DownCast(ent);
    if (!edgeList.IsNull()) {
      const Handle(TransferBRep_ShapeListBinder)& binder =
            Handle(TransferBRep_ShapeListBinder)::DownCast(TP->Find(ent));
      if (binder.IsNull()) continue;
      if (binder->NbShapes() != edgeList->NbEdges()) continue;

      for (Standard_Integer iedge = 1; iedge <= edgeList->NbEdges(); iedge++) {
        Handle(IGESData_IGESEntity) Curve = edgeList->Curve(iedge);
        S = binder->Edge(iedge);
        if (S.IsNull()) continue;
        EG_igesEntityColor(Curve, match, TP, S, shapeData);
      }
      continue;
    }

    // check each curve for an IGESGeom_CompositeCurve
    const Handle(IGESGeom_CompositeCurve)& compCurve =
          Handle(IGESGeom_CompositeCurve)::DownCast(ent);
    if (!compCurve.IsNull()) {
      const Handle(TransferBRep_ShapeBinder)& binder =
            Handle(TransferBRep_ShapeBinder)::DownCast(TP->Find(ent));
      if (binder.IsNull()) continue;

      TopoDS_Shape wire = binder->Result();
      if (wire.IsNull()) continue;
      match.Find(wire);

      TopExp_Explorer Exp;
      Exp.Init(wire, TopAbs_EDGE);

      for (Standard_Integer iedge = 1;
          iedge <= compCurve->NbCurves(); iedge++, Exp.Next()) {
        Handle(IGESData_IGESEntity) Curve = compCurve->Curve(iedge);
        S = Exp.Current();
        EG_igesEntityColor(Curve, match, TP, S, shapeData);
      }
      continue;
    }
  }
}


static void
EG_importScale(const char *reader, const char *units, double *scale,
               const char **wunits)
{
  if ((strcasecmp(units,"inch")   == 0) ||
      (strcasecmp(units,"inches") == 0) ||
      (strcasecmp(units,"in")     == 0)) {
    *scale = 1.0/25.4;
    if (wunits != NULL) *wunits = "INCH";
    printf("  %s Info: Using %s\n", reader, units);
  } else if ((strcasecmp(units,"foot") == 0) ||
             (strcasecmp(units,"feet") == 0) ||
             (strcasecmp(units,"ft")   == 0)) {
    *scale = 1.0/304.8;
    if (wunits != NULL) *wunits = "FT";
    printf("  %s Info: Using %s\n", reader, units);
  } else if ((strcasecmp(units,"metre")  == 0) ||
             (strcasecmp(units,"meter")  == 0) ||
             (strcasecmp(units,"meters") == 0) ||
             (strcasecmp(units,"m")      == 0)) {
    *scale = 1.0/1000.0;
    if (wunits != NULL) *wunits = "M";
    printf("  %s Info: Using %s\n", reader, units);
  } else if ((strcasecmp(units,"centimetre")  == 0) ||
             (strcasecmp(units,"centimeter")  == 0) ||
             (strcasecmp(units,"centimeters") == 0) ||
             (strcasecmp(units,"cm")          == 0)) {
    *scale = 1.0/10.0;
    if (wunits != NULL) *wunits = "CM";
    printf("  %s Info: Using %s\n", reader, units);
  } else if ((strcasecmp(units,"millimetre")  == 0) ||
             (strcasecmp(units,"millimeter")  == 0) ||
             (strcasecmp(units,"millimeters") == 0) ||
             (strcasecmp(units,"mm")          == 0)) {
    if (wunits != NULL) *wunits = "MM";
    printf("  %s Info: Using %s\n", reader, units);
  } else {
    if (wunits != NULL) *wunits = "MM";
    printf(" EGADS %s Info: Cannot convert %s -- using millimeters!\n",
           reader, units);
  }
}


int
EG_loadModel(egObject *context, int bflg, const char *name, egObject **model)
{
  int          i, j, stat, outLevel, len, nattr, nerr, hite, hitf, egads = 0;
  int          oclass, ibody, *invalid = NULL, nbs = 0;
  double       scale  = 1.0;
  egObject     *omodel, *aobj;
  TopoDS_Shape source;
  egadsModel   *mshape = NULL;
  egadsShapeData shapeData;
  Handle(TCollection_HAsciiString) units;
  FILE         *fp;

  *model = NULL;
  if (context == NULL)               return EGADS_NULLOBJ;
  if (context->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (context->oclass != CONTXT)     return EGADS_NOTCNTX;
  if (EG_sameThread(context))        return EGADS_CNTXTHRD;
  outLevel = EG_outLevel(context);

  if (name == NULL) {
    if (outLevel > 0)
      printf(" EGADS Warning: NULL Filename (EG_loadModel)!\n");
    return EGADS_NONAME;
  }

  /* does file exist? */

  fp = fopen(name, "r");
  if (fp == NULL) {
    if (outLevel > 0)
      printf(" EGADS Warning: File %s Not Found (EG_loadModel)!\n", name);
    return EGADS_NOTFOUND;
  }
  fclose(fp);

  /* find extension */

  len = strlen(name);
  for (i = len-1; i > 0; i--)
    if (name[i] == '.') break;
  if (i == 0) {
    if (outLevel > 0)
      printf(" EGADS Warning: No Extension in %s (EG_loadModel)!\n", name);
    return EGADS_NODATA;
  }

  if ((strcasecmp(&name[i],".step") == 0) ||
      (strcasecmp(&name[i],".stp") == 0)) {

    /* STEP files */

    STEPControl_Reader aReader;
    IFSelect_ReturnStatus status = aReader.ReadFile(name);
    if (status != IFSelect_RetDone) {
      if (outLevel > 0)
        printf(" EGADS Error: STEP Read of %s = %d (EG_loadModel)!\n",
               name, status);
      return EGADS_NOLOAD;
    }

    // inspect the root transfers
    if (outLevel > 2)
      aReader.PrintCheckLoad(Standard_False, IFSelect_ItemsByEntity);

    if ((bflg&16) != 0) egads = -1;
    if ((bflg&4)  == 0) {
      TColStd_SequenceOfAsciiString unitLength, unitAngle, solidAngle;
      aReader.FileUnits(unitLength, unitAngle, solidAngle);
      if (unitLength.Length() >= 1) {
        if (unitLength.Length() > 1)
          printf(" EGADS Info: # unitLengths = %d\n", unitLength.Length());
        units = new TCollection_HAsciiString(unitLength(1).ToCString());
        EG_importScale("STEP Reader", units->ToCString(), &scale, NULL);
      }
    }

    int nroot = aReader.NbRootsForTransfer();
    if (outLevel > 1)
      printf(" EGADS Info: %s Entries = %d\n", name, nroot);

    for (i = 1; i <= nroot; i++) {
      Standard_Boolean ok = aReader.TransferRoot(i);
      if ((!ok) && (outLevel > 0))
        printf(" EGADS Warning: Transfer %d/%d is not OK!\n", i, nroot);
    }

    nbs = aReader.NbShapes();
    if (nbs <= 0) {
      if (outLevel > 0)
        printf(" EGADS Error: %s has No Shapes (EG_loadModel)!\n",
               name);
      return EGADS_NOLOAD;
    }
    if (outLevel > 1)
      printf(" EGADS Info: %s has %d Shape(s)\n", name, nbs);

    const Handle(XSControl_WorkSession)& workSession = aReader.WS();
    const Handle(XSControl_TransferReader)& TR = workSession->TransferReader();
#ifdef STEPASSATTRS
    Handle(Standard_Type) tNAUO = STANDARD_TYPE(
                                       StepRepr_NextAssemblyUsageOccurrence);
    Handle(Standard_Type) tPD   = STANDARD_TYPE(StepBasic_ProductDefinition);
#endif

    // Get all the labels from the STEP file
    EG_labelsFromTransferReader(aReader.OneShape(), TR, shapeData);

    // Get colors from step file
    EG_stepColors(aReader.OneShape(), workSession, shapeData);

    TopoDS_Compound compound;
    BRep_Builder    builder3D;
    builder3D.MakeCompound(compound);
    for (i = 1; i <= nbs; i++) {
      TopoDS_Shape aShape = aReader.Shape(i);
#ifdef OLD_AND_SLOW_LABEL_FINDER
      TopTools_IndexedMapOfShape MapAll;
      TopExp::MapShapes(aShape, MapAll);

      for (int j = 1; j <= MapAll.Extent(); j++) {
        TopoDS_Shape aShape = MapAll(j);
        Handle(Standard_Transient) ent = TR->EntityFromShapeResult(aShape, 1);
        if (ent.IsNull())          ent = TR->EntityFromShapeResult(aShape,-1);
        if (ent.IsNull())          ent = TR->EntityFromShapeResult(aShape, 4);
        if (ent.IsNull())          continue;
        Handle(StepRepr_RepresentationItem) aReprItem;
        aReprItem = Handle(StepRepr_RepresentationItem)::DownCast(ent);
        if (aReprItem.IsNull())    continue;
        const char *STEPname = aReprItem->Name()->ToCString();
        if (STEPname == NULL)      continue;
        if (strlen(STEPname) == 0) continue;
        shapeData.Add(aShape, STEPname);
        /* printf(" %d/%d %s -> Name found: %s\n", j, MapAll.Extent(),
               TopAbs::ShapeTypeToString(aShape.ShapeType()), STEPname ); */
      }
#endif

#ifdef STEPASSATTRS
      Handle(Standard_Transient) ent = TR->EntityFromShapeResult(aShape, 3);
      if (!ent.IsNull()) {
        if (ent->DynamicType() == tNAUO) {
          Handle(StepRepr_NextAssemblyUsageOccurrence) NAUO =
              Handle(StepRepr_NextAssemblyUsageOccurrence)::DownCast(ent);
          if (!NAUO.IsNull()) {
            Interface_EntityIterator subs = workSession->Graph().Sharings(NAUO);
            for (subs.Start(); subs.More(); subs.Next()) {
              Handle(StepRepr_ProductDefinitionShape) PDS =
                Handle(StepRepr_ProductDefinitionShape)::DownCast(subs.Value());
              if (!PDS.IsNull()) {
                Handle(StepBasic_ProductDefinitionRelationship) PDR =
                              PDS->Definition().ProductDefinitionRelationship();
                if (!PDR.IsNull()) {
                  if (PDR->HasDescription() && PDR->Description()->Length() > 0) {
                    printf(" NAUOa = %s\n", PDR->Description()->ToCString());
                  } else if (PDR->Name()->Length() > 0) {
                    printf(" NAUOb = %s\n", PDR->Name()->ToCString());
                  } else {
                    printf(" NAUOc = %s\n", PDR->Id()->ToCString());
                  }
                }
              }
            }
          }
        } else if (ent->DynamicType() == tPD) {
          Handle(StepBasic_ProductDefinition) PD =
              Handle(StepBasic_ProductDefinition)::DownCast(ent);
          if (!PD.IsNull()) {
            Handle(StepBasic_Product) Prod = PD->Formation()->OfProduct();
            if (Prod->Name()->UsefullLength() > 0) {
              printf(" PDa   = %s\n", Prod->Name()->ToCString());
            } else {
              printf(" PDb   = %s\n", Prod->Id()->ToCString());
            }
          }
        }
      }
#endif
      if ((bflg&8) != 0) {
        ShapeUpgrade_UnifySameDomain uShape(aShape, Standard_True,
                                            Standard_True, Standard_True);
//      uShape.SetLinearTolerance(100.0*Precision::Confusion());
//      uShape.SetAngularTolerance(10.0*Precision::Angular());
        uShape.Build();
        aShape = shapeData.Update(aShape, uShape.Shape(), uShape.History());
      }
      if (scale != 1.0) {
        gp_Trsf form = gp_Trsf();
        form.SetValues(scale, 0.0,   0.0,   0.0,
                       0.0,   scale, 0.0,   0.0,
                       0.0,   0.0,   scale, 0.0);
        BRepBuilderAPI_Transform xForm(aShape, form, Standard_True);
        if (!xForm.IsDone()) {
          printf(" EGADS Warning: Can't scale Body %d (EG_loadModel)!\n", i);
        } else {
          aShape = shapeData.Update(aShape, xForm);
        }
      }

      if (aShape.ShapeType() == TopAbs_COMPOUND) {
        // Wires come in as composite of edges...
        TopExp_Explorer Exp;
        BRepBuilderAPI_MakeWire MW;
        for (Exp.Init(aShape, TopAbs_EDGE,  TopAbs_WIRE);
             Exp.More(); Exp.Next()) {
          MW.Add( TopoDS::Edge(Exp.Current()));
        }
        if (MW.IsDone()) {
          TopoDS_Shape Wire = shapeData.Update(aShape, MW.Wire());
          TopExp_Explorer WExp(Wire, TopAbs_EDGE);
          for (Exp.Init(aShape, TopAbs_EDGE,  TopAbs_WIRE);
               Exp.More() && WExp.More(); Exp.Next(), WExp.Next()) {
            shapeData.Update(Exp.Current(), WExp.Current());
          }
          builder3D.Add(compound, Wire);
        }
        for (Exp.Init(aShape, TopAbs_WIRE,  TopAbs_FACE);
             Exp.More(); Exp.Next()) builder3D.Add(compound, Exp.Current());
        for (Exp.Init(aShape, TopAbs_FACE,  TopAbs_SHELL);
             Exp.More(); Exp.Next()) builder3D.Add(compound, Exp.Current());
        for (Exp.Init(aShape, TopAbs_SHELL, TopAbs_SOLID);
             Exp.More(); Exp.Next()) builder3D.Add(compound, Exp.Current());
        for (Exp.Init(aShape, TopAbs_SOLID);
             Exp.More(); Exp.Next()) builder3D.Add(compound, Exp.Current());
      } else {
        builder3D.Add(compound, aShape);
      }
    }
    source = compound;

    if (outLevel > 0) {
      hite = hitf = 0;
      TopTools_IndexedMapOfShape MapE, MapF;
      TopExp::MapShapes(source, TopAbs_EDGE, MapE);
      for (i = 1; i <= MapE.Extent(); i++) {
        Standard_Real t1, t2;
        TopoDS_Shape shape = MapE(i);
        TopoDS_Edge  Edge  = TopoDS::Edge(shape);
        Handle(Geom_Curve) hCurve = BRep_Tool::Curve(Edge, t1, t2);
        Handle(Geom_BSplineCurve) hBSpline =
          Handle(Geom_BSplineCurve)::DownCast(hCurve);
        if (hBSpline.IsNull()) continue;
        if (hBSpline->Continuity() == GeomAbs_C0)  hite++;
      }
      TopExp::MapShapes(source, TopAbs_FACE, MapF);
      for (int i = 1; i <= MapF.Extent(); i++) {
        TopoDS_Shape shape = MapF(i);
        TopoDS_Face  Face  = TopoDS::Face(shape);
        Handle(Geom_Surface) hSurface = BRep_Tool::Surface(Face);
        Handle(Geom_BSplineSurface) hBSpline =
          Handle(Geom_BSplineSurface)::DownCast(hSurface);
        if (hBSpline.IsNull()) continue;
        if (hBSpline->Continuity() == GeomAbs_C0)  hitf++;
      }
      if (hite+hitf != 0)
        printf(" EGADS Info: Import Model has %d C0 Edges & %d C0 Faces!\n",
               hite, hitf);
    }

  } else if ((strcasecmp(&name[i],".iges") == 0) ||
             (strcasecmp(&name[i],".igs") == 0)) {

    /* IGES files */

    IGESControl_Reader iReader;
    Standard_Integer stats = iReader.ReadFile(name);
    if (stats != IFSelect_RetDone) {
      if (outLevel > 0)
        printf(" EGADS Error: IGES Read of %s = %d (EG_loadModel)!\n",
               name, stats);
      return EGADS_NOLOAD;
    }
    if ((bflg&4) == 0) {
      Handle(IGESData_IGESModel) aModel =
                         Handle(IGESData_IGESModel)::DownCast(iReader.Model());
      if(!aModel.IsNull()) {
        units = aModel->GlobalSection().UnitName();
        if (!units.IsNull())
          EG_importScale("IGES Reader", units->ToCString(), &scale, NULL);
      }
    }
    iReader.TransferRoots();
    if ((bflg&16) != 0) egads = -1;

    nbs = iReader.NbShapes();
    if (nbs <= 0) {
      if (outLevel > 0)
        printf(" EGADS Error: %s has No Shapes (EG_loadModel)!\n",
               name);
      return EGADS_NOLOAD;
    }
    if (outLevel > 1)
      printf(" EGADS Info: %s has %d Shape(s)\n", name, nbs);

    const Handle(XSControl_WorkSession)& workSession = iReader.WS();
    const Handle(XSControl_TransferReader)& transferReader =
          workSession->TransferReader();

    // Get all the labels from the IGES file
    EG_labelsFromTransferReader(iReader.OneShape(), transferReader, shapeData);

    // Get all the colors from the IGES file
    EG_igesColors(iReader.OneShape(), workSession, shapeData);

    TopoDS_Compound compound;
    BRep_Builder    builder3D;
    builder3D.MakeCompound(compound);
    for (i = 1; i <= nbs; i++) {
      TopoDS_Shape aShape = iReader.Shape(i);

#ifdef OLD_AND_SLOW_LABEL_FINDER
      const Handle(IGESData_IGESEntity)& shapeEntity =
            Handle(IGESData_IGESEntity)::DownCast(
                  transferReader->EntityFromShapeResult(aShape, -1));
      if (shapeEntity->HasName())
        shapeData.Add(aShape, shapeEntity->NameValue()->ToCString());
#endif

      if ((bflg&8) != 0) {
        ShapeUpgrade_UnifySameDomain uShape(aShape);
        uShape.Build();
        aShape = uShape.Shape();
      }
      if (scale != 1.0) {
        gp_Trsf form = gp_Trsf();
        form.SetValues(scale, 0.0,   0.0,   0.0,
                       0.0,   scale, 0.0,   0.0,
                       0.0,   0.0,   scale, 0.0);
        BRepBuilderAPI_Transform xForm(aShape, form, Standard_True);
        if (!xForm.IsDone()) {
          printf(" EGADS Warning: Can't scale Body %d (EG_loadModel)!\n", i);
        } else {
          aShape = shapeData.Update(aShape, xForm.ModifiedShape(aShape));
        }
      }
      builder3D.Add(compound, aShape);
    }
    source = compound;

    if (outLevel > 0) {
      hite = hitf = 0;
      TopTools_IndexedMapOfShape MapE, MapF;
      TopExp::MapShapes(source, TopAbs_EDGE, MapE);
      for (i = 1; i <= MapE.Extent(); i++) {
        Standard_Real t1, t2;
        TopoDS_Shape shape = MapE(i);
        TopoDS_Edge  Edge  = TopoDS::Edge(shape);
        Handle(Geom_Curve) hCurve = BRep_Tool::Curve(Edge, t1, t2);
        Handle(Geom_BSplineCurve) hBSpline =
          Handle(Geom_BSplineCurve)::DownCast(hCurve);
        if (hBSpline.IsNull()) continue;
        if (hBSpline->Continuity() == GeomAbs_C0)  hite++;
      }
      TopExp::MapShapes(source, TopAbs_FACE, MapF);
      for (int i = 1; i <= MapF.Extent(); i++) {
        TopoDS_Shape shape = MapF(i);
        TopoDS_Face  Face  = TopoDS::Face(shape);
        Handle(Geom_Surface) hSurface = BRep_Tool::Surface(Face);
        Handle(Geom_BSplineSurface) hBSpline =
          Handle(Geom_BSplineSurface)::DownCast(hSurface);
        if (hBSpline.IsNull()) continue;
        if (hBSpline->Continuity() == GeomAbs_C0)  hitf++;
      }
      if (hite+hitf != 0)
        printf(" EGADS Info: Import Model has %d C0 Edges & %d C0 Faces!\n",
               hite, hitf);
    }

  } else if ((strcasecmp(&name[i],".brep") == 0) ||
             (strcasecmp(&name[i],".egads") == 0)) {

    /* Native OCC file */
    if (strcasecmp(&name[i],".egads") == 0) egads = 1;

    BRep_Builder builder;
    if (!BRepTools::Read(source, name, builder)) {
      if (outLevel > 0)
        printf(" EGADS Warning: Read Error on %s (EG_loadModel)!\n", name);
      return EGADS_NOLOAD;
    }

  } else {
    if (outLevel > 0)
      printf(" EGADS Warning: Extension in %s Not Supported (EG_loadModel)!\n",
             name);
    return EGADS_NODATA;
  }

  int nWire  = 0;
  int nFace  = 0;
  int nSheet = 0;
  int nSolid = 0;

  TopExp_Explorer Exp;
  for (Exp.Init(source, TopAbs_WIRE,  TopAbs_FACE);
       Exp.More(); Exp.Next()) nWire++;
  for (Exp.Init(source, TopAbs_FACE,  TopAbs_SHELL);
       Exp.More(); Exp.Next()) nFace++;
  for (Exp.Init(source, TopAbs_SHELL, TopAbs_SOLID);
       Exp.More(); Exp.Next()) nSheet++;
  for (Exp.Init(source, TopAbs_SOLID); Exp.More(); Exp.Next()) nSolid++;

  if (outLevel > 1)
    printf("\n EGADS Info: %s has %d Solids, %d Sheets, %d Faces and %d Wires\n",
           name, nSolid, nSheet, nFace, nWire);

  int nBody = nWire+nFace+nSheet+nSolid;
  /* can we promote Edges to Wires? */
  if ((nBody == 0) || (egads == -1)) {
    double tol = 1.e-7;
/*
    for (Exp.Init(source, TopAbs_VERTEX); Exp.More(); Exp.Next()) {
      TopoDS_Vertex Vert = TopoDS::Vertex(Exp.Current());
      double tolv = BRep_Tool::Tolerance(Vert);
      if (tolv > tol) tol = tolv;
    }
 */
    int nEdge, add = 0;
    j = 0;
    for (Exp.Init(source, TopAbs_EDGE, TopAbs_WIRE); Exp.More(); Exp.Next()) j++;
    if (j != 0) {
      nEdge = j;
      TopoDS_Edge *Edges = new TopoDS_Edge[nEdge];
      /* remove small Edges */
      j = 0;
      for (Exp.Init(source, TopAbs_EDGE, TopAbs_WIRE); Exp.More(); Exp.Next()) {
        TopoDS_Shape shape = Exp.Current();
        TopoDS_Edge  Edge  = TopoDS::Edge(shape);
        Standard_Real t1, t2;
        Handle(Geom_Curve) hCurve = BRep_Tool::Curve(Edge, t1, t2);
        GeomAdaptor_Curve AC(hCurve);
        if (GCPnts_AbscissaPoint::Length(AC, t1, t2) < tol) continue;
        Edges[j] = Edge;
        j++;
      }
      nEdge = j;
      TopoDS_Compound compound;
      BRep_Builder    builder3D;
      builder3D.MakeCompound(compound);
      if (nBody != 0) builder3D.Add(compound, source);
      while (j != 0) {
        for (i = 0; i < nEdge; i++)
          if (!Edges[i].IsNull()) {
            BRepBuilderAPI_MakeWire MW;
            try {
              MW.Add(Edges[i]);
              TopoDS_Vertex V1, V2;
              TopExp::Vertices(Edges[i], V2, V1, Standard_True);
              if (outLevel > 1) printf(" start Edge %d:", i+1);
              Edges[i].Nullify();
              if (V2.IsSame(V1)) {
                if (outLevel > 1) printf("\n");
                break;
              }
              gp_Pnt pv1 = BRep_Tool::Pnt(V1);
              gp_Pnt pv2 = BRep_Tool::Pnt(V2);
              int hit;
              do {
                hit = 0;
                for (int k = 0; k < nEdge; k++)
                  if (!Edges[k].IsNull()) {
                    TopExp::Vertices(Edges[k], V2, V1, Standard_True);
                    gp_Pnt pv = BRep_Tool::Pnt(V1);
                    if (pv.Distance(pv1) < tol) {
                      if (outLevel > 1) printf(" --%d", k+1);
                      MW.Add(Edges[k]);
                      Edges[k].Nullify();
                      pv1 = pv;
                      hit++;
                      break;
                    }
                    if (pv.Distance(pv2) < tol) {
                      if (outLevel > 1) printf(" -+%d", k+1);
                      MW.Add(Edges[k]);
                      Edges[k].Nullify();
                      pv2 = pv;
                      hit++;
                      break;
                    }
                    pv = BRep_Tool::Pnt(V2);
                    if (pv.Distance(pv1) < tol) {
                      if (outLevel > 1) printf(" +-%d", k+1);
                      MW.Add(Edges[k]);
                      Edges[k].Nullify();
                      pv1 = pv;
                      hit++;
                      break;
                    }
                    if (pv.Distance(pv2) < tol) {
                      if (outLevel > 1) printf(" ++%d", k+1);
                      MW.Add(Edges[k]);
                      Edges[k].Nullify();
                      pv2 = pv;
                      hit++;
                      break;
                    }
                  }
              } while (hit != 0);
              if (outLevel > 1) printf("\n");
            }
            catch (const Standard_Failure& e) {
              printf(" EGADS Warning: Cannot Add Edge in Wire!\n");
              printf("                %s\n", e.GetMessageString());
              continue;
            }
            catch (...) {
              printf(" EGADS Warning: Cannot Add Edge in Wire!\n");
              continue;
            }
            if (MW.Error()) {
              if (outLevel > 0)
                printf(" EGADS Error: Problem with Imported Edge!\n");
              continue;
            }
            if (!MW.IsDone()) {
              if (outLevel > 0)
                printf(" EGADS Error: Problem with Loop Conversion!\n");
              continue;
            }
            TopoDS_Wire wire = MW.Wire();
            builder3D.Add(compound, wire);
            nWire++;
            nBody++;
            add++;
          }
        j = 0;
        for (i = 0; i < nEdge; i++) if (!Edges[i].IsNull()) j++;
      }
      source = compound;
      delete [] Edges;
    }
    if (add != 0)
      if (outLevel > 0)
        printf(" EGADS Info: %d Edges converted to %d WireBodies (EG_loadModel)!\n",
               nEdge, add);
  }
  if (nBody == 0) {
    source.Nullify();
    if (outLevel > 0)
      printf(" EGADS Warning: Nothing found in %s (EG_loadModel)!\n", name);
    return EGADS_NODATA;
  }
  if (egads == 0) {
    invalid = (int *) EG_alloc(nBody*sizeof(int));
    if (invalid == NULL) {
      printf(" EGADS Warning: Cannot check for Invalid Bodies (EG_loadModel)!\n");
    } else {
      for (i = 0; i < nBody; i++) invalid[i] = 0;
    }
  }

  mshape              = new egadsModel;
  mshape->shape       = source;
  mshape->nobjs       = nBody;
  mshape->nbody       = nBody;
  mshape->bbox.filled = 0;
  mshape->bodies      = new egObject*[nBody];
  for (i = 0; i < nBody; i++) {
    stat = EG_makeObject(context, &mshape->bodies[i]);
    if (stat != EGADS_SUCCESS) {
      for (int j = 0; j < i; j++) {
        egObject  *obj   = mshape->bodies[j];
        egadsBody *pbody = (egadsBody *) obj->blind;
        delete pbody;
        EG_deleteObject(mshape->bodies[j]);
      }
      delete [] mshape->bodies;
      delete mshape;
      return stat;
    }
    egObject  *pobj    = mshape->bodies[i];
    egadsBody *pbody   = new egadsBody;
    pbody->nodes.objs  = NULL;
    pbody->edges.objs  = NULL;
    pbody->loops.objs  = NULL;
    pbody->faces.objs  = NULL;
    pbody->shells.objs = NULL;
    pbody->senses      = NULL;
    pbody->bbox.filled = 0;
    pbody->massFill    = 0;
    pobj->blind        = pbody;
  }

  i = 0;
  for (Exp.Init(mshape->shape, TopAbs_WIRE,  TopAbs_FACE);
       Exp.More(); Exp.Next()) {
    egObject  *obj   = mshape->bodies[i++];
    egadsBody *pbody = (egadsBody *) obj->blind;
    pbody->shape     = Exp.Current();
    if (egads == 0) {
      BRepCheck_Analyzer wCheck(pbody->shape);
      if (!wCheck.IsValid()) {
        Handle_ShapeFix_Shape sfs = new ShapeFix_Shape(pbody->shape);
        sfs->Perform();
        TopoDS_Shape fixedShape = sfs->Shape();
        if (fixedShape.IsNull()) {
          if (outLevel > 0)
            printf(" EGADS Warning: Body %d is an Invalid WireBody!\n", i);
          if (invalid != NULL) invalid[i-1] = 1;
        } else {
          BRepCheck_Analyzer wfCheck(fixedShape);
          if (!wfCheck.IsValid()) {
            if (outLevel > 0)
              printf(" EGADS Warning: Fixed Body %d is an Invalid WireBody!\n",
                     i);
            if (invalid != NULL) invalid[i-1] = 1;
          } else {
            pbody->shape = shapeData.Update(pbody->shape, fixedShape);
          }
        }
      }
    }
  }
  for (Exp.Init(mshape->shape, TopAbs_FACE,  TopAbs_SHELL);
       Exp.More(); Exp.Next()) {
    egObject  *obj   = mshape->bodies[i++];
    egadsBody *pbody = (egadsBody *) obj->blind;
    pbody->shape     = Exp.Current();
    if (egads == 0) {
      BRepCheck_Analyzer fCheck(pbody->shape);
      if (!fCheck.IsValid()) {
        Handle_ShapeFix_Shape sfs = new ShapeFix_Shape(pbody->shape);
        sfs->Perform();
        TopoDS_Shape fixedShape = sfs->Shape();
        if (fixedShape.IsNull()) {
          if (outLevel > 0)
            printf(" EGADS Warning: Body %d is an Invalid FaceBody!\n", i);
          if (invalid != NULL) invalid[i-1] = 1;
        } else {
          BRepCheck_Analyzer ffCheck(fixedShape);
          if (!ffCheck.IsValid()) {
            if (outLevel > 0)
              printf(" EGADS Warning: Fixed Body %d is an Invalid FaceBody!\n",
                     i);
            if (invalid != NULL) invalid[i-1] = 1;
          } else {
            pbody->shape = shapeData.Update(pbody->shape, fixedShape);
          }
        }
      }
    }
  }
  for (Exp.Init(mshape->shape, TopAbs_SHELL, TopAbs_SOLID);
       Exp.More(); Exp.Next()) {
    egObject  *obj   = mshape->bodies[i++];
    egadsBody *pbody = (egadsBody *) obj->blind;
    pbody->shape     = Exp.Current();
    if (egads == 0) {
      BRepCheck_Analyzer sCheck(pbody->shape);
      if (!sCheck.IsValid()) {
        Handle_ShapeFix_Shape sfs = new ShapeFix_Shape(pbody->shape);
        sfs->Perform();
        TopoDS_Shape fixedShape = sfs->Shape();
        if (fixedShape.IsNull()) {
          if (outLevel > 0)
            printf(" EGADS Warning: Body %d is an Invalid SheetBody!\n", i);
          if (invalid != NULL) invalid[i-1] = 1;
        } else {
          if (fixedShape.ShapeType() != TopAbs_SHELL) {
            if (outLevel > 0)
              printf(" EGADS Reject: Fixed Body %d is No longer a SheetBody!\n",
                     i);
            if (invalid != NULL) invalid[i-1] = 1;
          } else {
            BRepCheck_Analyzer sfCheck(fixedShape);
            if (!sfCheck.IsValid()) {
              if (outLevel > 0)
                printf(" EGADS Warning: Fixed Body %d is an Invalid SheetBody!\n",
                       i);
              if (invalid != NULL) invalid[i-1] = 1;
            } else {
              pbody->shape = shapeData.Update(pbody->shape, fixedShape);
            }
          }
        }
      }
    }
  }
  for (Exp.Init(mshape->shape, TopAbs_SOLID); Exp.More(); Exp.Next()) {
    egObject  *obj   = mshape->bodies[i++];
    egadsBody *pbody = (egadsBody *) obj->blind;
    pbody->shape     = Exp.Current();
    if (egads == 0) {
      BRepCheck_Analyzer sCheck(pbody->shape);
      if (!sCheck.IsValid()) {
        Handle_ShapeFix_Shape sfs = new ShapeFix_Shape(pbody->shape);
        sfs->Perform();
        TopoDS_Shape fixedShape = sfs->Shape();
        if (fixedShape.IsNull()) {
          if (outLevel > 0)
            printf(" EGADS Warning: Body %d is an Invalid SolidBody!\n", i);
          if (invalid != NULL) invalid[i-1] = 1;
        } else {
          if (fixedShape.ShapeType() != TopAbs_SOLID) {
            if (outLevel > 0)
              printf(" EGADS Reject: Fixed Body %d is No longer a SolidBody!\n",
                     i);
            if (invalid != NULL) invalid[i-1] = 1;
          } else {
            BRepCheck_Analyzer sfCheck(fixedShape);
            if (!sfCheck.IsValid()) {
              if (outLevel > 0)
                printf(" EGADS Warning: Fixed Body %d is an Invalid SolidBody!\n",
                       i);
              if (invalid != NULL) invalid[i-1] = 1;
            } else {
              pbody->shape = shapeData.Update(pbody->shape, fixedShape);
            }
          }
        }
      }
    }
  }

  stat = EG_makeObject(context, &omodel);
  if (stat != EGADS_SUCCESS) {
    source.Nullify();
    for (i = 0; i < nBody; i++) {
      egObject  *obj   = mshape->bodies[i];
      egadsBody *pbody = (egadsBody *) obj->blind;
      delete pbody;
      EG_deleteObject(mshape->bodies[i]);
    }
    delete [] mshape->bodies;
    delete mshape;
    if (invalid != NULL) EG_free(invalid);
    return stat;
  }
  omodel->oclass = MODEL;
  omodel->blind  = mshape;
  EG_referenceObject(omodel, context);

  for (j = i = 0; i < nBody; i++) {
    egObject  *pobj  = mshape->bodies[i];
    egadsBody *pbody = (egadsBody *) pobj->blind;
    pobj->topObj     = omodel;
    if (((bflg&1) == 0) && (egads == 0)) {
      int ibi = i+1;
      if ((bflg&32) != 0) ibi = -ibi;           // undocumented -- force split
      EG_splitPeriodics(ibi, pbody, shapeData);
    }
    if (((bflg&2) != 0) && (egads == 0))
      EG_splitMultiplicity(pbody, shapeData, outLevel);
    stat = EG_traverseBody(context, i, pobj, omodel, pbody, &nerr);
    if (stat != EGADS_SUCCESS) {
      if (egads == 1) {
        mshape->nbody = i;
        EG_destroyTopology(omodel);
        if (invalid != NULL) EG_free(invalid);
        return stat;
      }
      if (outLevel > 0)
        printf(" EGADS Warning: Body  %d parse = %d!\n", i+1, stat);
      continue;
    }
    if (invalid != NULL)
      if (invalid[i] == 1)
        EG_attributeAdd(pobj, ".invalid", ATTRSTRING, 1, NULL, NULL, name);
    mshape->bodies[j] = mshape->bodies[i];
    j++;
  }
  if (j != mshape->nbody) {
    if (outLevel > 0)
      printf(" EGADS Inform:  %d Bodies not included!\n\n", mshape->nbody-j);
    mshape->nbody = j;
    if (j == 0) {
      EG_destroyTopology(omodel);
      if (invalid != NULL) EG_free(invalid);
      return EGADS_TOPOERR;
    }
  }
  *model = omodel;

  /* possibly assign units when doing an IGES/STEP read */
  if (!units.IsNull()) {
    for (int ibody = 0; ibody < mshape->nbody; ibody++) {
      egObject *pobj = mshape->bodies[ibody];
      EG_attributeAdd(pobj, ".lengthUnits", ATTRSTRING, 1, NULL, NULL,
                      units->ToCString());
    }
  }

  /* possibly assign attributes from IGES/STEP read */
  if (shapeData.LabelExtent() > 0) {
    for (int ibody = 0; ibody < mshape->nbody; ibody++) {
      egObject  *pobj   = mshape->bodies[ibody];
      egadsBody *pbody  = (egadsBody *) pobj->blind;
      i = shapeData.LabelFindIndex(pbody->shape);
      if (i > 0) {
        const char *value = shapeData.label(i).shapeName;
        EG_attributeAdd(pobj, "Name", ATTRSTRING, 1, NULL, NULL, value);
      }
      shapeData.AddShapeName(pbody->nodes);
      shapeData.AddShapeName(pbody->edges);
      shapeData.AddShapeName(pbody->loops);
      shapeData.AddShapeName(pbody->faces);
      shapeData.AddShapeName(pbody->shells);
    }
  }
  if (shapeData.BodyColorExtent() > 0) {
    TopTools_IndexedMapOfShape bodyMap;
    for (int ibody = 0; ibody < mshape->nbody; ibody++) {
      egObject  *pobj   = mshape->bodies[ibody];
      egadsBody *pbody  = (egadsBody *) pobj->blind;
      bodyMap.Add(pbody->shape);
    }

    for (i = 1; i <= shapeData.BodyColorExtent(); i++) {
      const Quantity_Color* coloredge = shapeData.EdgeColor(i);
      const Quantity_Color* colorface = shapeData.FaceColor(i);
      const TopoDS_Shape& aShape = shapeData.BodyColorFindKey(i);

      ibody = bodyMap.FindIndex(aShape);
      if (ibody == 0) continue;
      egObject *pobj = mshape->bodies[ibody-1];

      if (coloredge != NULL) {
        const double rgb[3] = {coloredge->Red(), coloredge->Green(), coloredge->Blue()};

        const char* name = NULL;
        if (rgb[0] == 1   && rgb[1] == 0   && rgb[2] == 0  ) name = "red";
        else if (rgb[0] == 0   && rgb[1] == 1   && rgb[2] == 0  ) name = "green";
        else if (rgb[0] == 0   && rgb[1] == 0   && rgb[2] == 1  ) name = "blue";
        else if (rgb[0] == 1   && rgb[1] == 1   && rgb[2] == 0  ) name = "yellow";
        else if (rgb[0] == 1   && rgb[1] == 0   && rgb[2] == 1  ) name = "magenta";
        else if (rgb[0] == 0   && rgb[1] == 1   && rgb[2] == 1  ) name = "cyan";
        else if (rgb[0] == 1   && rgb[1] == 1   && rgb[2] == 1  ) name = "white";
        else if (rgb[0] == 0   && rgb[1] == 0   && rgb[2] == 0  ) name = "black";
        else if (rgb[0] == 1   && rgb[1] == 0.5 && rgb[2] == 0.5) name = "lred";
        else if (rgb[0] == 0.5 && rgb[1] == 1   && rgb[2] == 0.5) name = "lgreen";
        else if (rgb[0] == 0.5 && rgb[1] == 0.5 && rgb[2] == 1  ) name = "lblue";

        if (name != NULL) {
          EG_attributeAdd(pobj, "ColorEdge", ATTRSTRING, 1, NULL, NULL, name);
        } else {
          EG_attributeAdd(pobj, "ColorEdge", ATTRREAL, 3, NULL, rgb, NULL);
        }
      }
      if (colorface != NULL) {
        const double rgb[3] = {colorface->Red(), colorface->Green(), colorface->Blue()};

        const char* name = NULL;
        if (rgb[0] == 1   && rgb[1] == 0   && rgb[2] == 0  ) name = "red";
        else if (rgb[0] == 0   && rgb[1] == 1   && rgb[2] == 0  ) name = "green";
        else if (rgb[0] == 0   && rgb[1] == 0   && rgb[2] == 1  ) name = "blue";
        else if (rgb[0] == 1   && rgb[1] == 1   && rgb[2] == 0  ) name = "yellow";
        else if (rgb[0] == 1   && rgb[1] == 0   && rgb[2] == 1  ) name = "magenta";
        else if (rgb[0] == 0   && rgb[1] == 1   && rgb[2] == 1  ) name = "cyan";
        else if (rgb[0] == 1   && rgb[1] == 1   && rgb[2] == 1  ) name = "white";
        else if (rgb[0] == 0   && rgb[1] == 0   && rgb[2] == 0  ) name = "black";
        else if (rgb[0] == 1   && rgb[1] == 0.5 && rgb[2] == 0.5) name = "lred";
        else if (rgb[0] == 0.5 && rgb[1] == 1   && rgb[2] == 0.5) name = "lgreen";
        else if (rgb[0] == 0.5 && rgb[1] == 0.5 && rgb[2] == 1  ) name = "lblue";

        if (name != NULL) {
          EG_attributeAdd(pobj, "ColorFace", ATTRSTRING, 1, NULL, NULL, name);
        } else {
          EG_attributeAdd(pobj, "ColorFace", ATTRREAL, 3, NULL, rgb, NULL);
        }
      }
    }
  }
  if (shapeData.ColorExtent() > 0) {
    for (int ibody = 0; ibody < mshape->nbody; ibody++) {
      egObject  *pobj   = mshape->bodies[ibody];
      egadsBody *pbody  = (egadsBody *) pobj->blind;

      shapeData.AddShapeColor(pbody->nodes);
      shapeData.AddShapeColor(pbody->edges);
      shapeData.AddShapeColor(pbody->loops);
      shapeData.AddShapeColor(pbody->faces);
      shapeData.AddShapeColor(pbody->shells);
    }
  }

  if (invalid != NULL) EG_free(invalid);
  if (egads != 1) return EGADS_SUCCESS;

  /* get the attributes from the EGADS files */

  fp = fopen(name, "r");
  if (fp == NULL) {
    printf(" EGADS Info: Cannot reOpen %s (EG_loadModel)!\n", name);
    return EGADS_SUCCESS;
  }
  char line[81];
  for (;;) {
    line[0] = line[1] = ' ';
    if (fgets(line, 81, fp) == NULL) break;
    if ((line[0] == '#') && (line[1] == '#')) break;
  }

  // got the header
  if ((line[0] == '#') && (line[1] == '#')) {
    if (outLevel > 1) printf(" Header = %s\n", line);
    // get number of model attributes
    fscanf(fp, "%d", &nattr);
    if (nattr != 0) EG_readAttrs(omodel, nattr, fp);
    for (i = 0; i < nBody; i++) {
      int otype,  oindex;
      int rsolid, rshell, rface, rloop, redge, rnode;
      int nsolid, nshell, nface, nloop, nedge, nnode;

      fscanf(fp, " %d %d %d %d %d %d %d", &rsolid, &rshell,
             &rface, &rloop, &redge, &rnode, &nattr);
      if (outLevel > 2)
        printf(" read = %d %d %d %d %d %d %d\n", rsolid, rshell,
               rface, rloop, redge, rnode, nattr);
      egObject  *pobj  = mshape->bodies[i];
      egadsBody *pbody = (egadsBody *) pobj->blind;
      nnode  = pbody->nodes.map.Extent();
      nedge  = pbody->edges.map.Extent();
      nloop  = pbody->loops.map.Extent();
      nface  = pbody->faces.map.Extent();
      nshell = pbody->shells.map.Extent();
      nsolid = 0;
      if (pobj->mtype == SOLIDBODY) nsolid = 1;
      if ((nnode != rnode) || (nedge  != redge)  || (nloop  != rloop) ||
          (nface != rface) || (nshell != rshell) || (nsolid != rsolid)) {
        printf(" EGADS Info: %d %d, %d %d, %d %d, %d %d, %d %d, %d %d",
               nnode, rnode, nedge,  redge,  nloop,  rloop,
               nface, rface, nshell, rshell, nsolid, rsolid);
        printf("  MisMatch on Attributes (EG_loadModel)!\n");
        fclose(fp);
        return EGADS_SUCCESS;
      }
      // got the correct body -- transfer the attributes
      if (nattr != 0) EG_readAttrs(pobj, nattr, fp);
      for (;;)  {
        j = fscanf(fp, "%d %d %d\n", &otype, &oindex, &nattr);
        if (outLevel > 2)
          printf(" %d:  attr header = %d %d %d\n",
                 j, otype, oindex, nattr);
        if (j     != 3) break;
        if (otype == 0) break;
        if (otype == 1) {
          aobj = pbody->shells.objs[oindex];
        } else if (otype == 2) {
          aobj = pbody->faces.objs[oindex];
        } else if (otype == 3) {
          aobj = pbody->loops.objs[oindex];
        } else if (otype == 4) {
          aobj = pbody->edges.objs[oindex];
        } else {
          aobj = pbody->nodes.objs[oindex];
        }
        EG_readAttrs(aobj, nattr, fp);
      }
    }

    /* get the ancillary objects from the EGADS files */
    line[0] = line[1] = ' ';
    fgets(line, 81, fp);
    if ((line[0] == '#') && (line[1] == '#')) {
      j = fscanf(fp, "%hd", &omodel->mtype);
      if ((j != 1) || (omodel->mtype < mshape->nbody)) {
        printf(" EGADS Info: Ext failure in %s  %d %d (EG_loadModel)!\n",
               name, j, omodel->mtype);
        omodel->mtype = 0;
        fclose(fp);
        return EGADS_SUCCESS;
      }
      egObject** bodies = new egObject*[omodel->mtype];
      for (j = 0; j < mshape->nbody; j++) bodies[j] = mshape->bodies[j];
      for (j = mshape->nbody; j < omodel->mtype; j++) bodies[j] = NULL;
      delete [] mshape->bodies;
      mshape->nobjs  = omodel->mtype;
      mshape->bodies = bodies;
      for (j = mshape->nbody; j < omodel->mtype; j++) {
        i = fscanf(fp, "%d %d", &oclass, &ibody);
        if (i != 2) {
          omodel->mtype = j;
          mshape->nobjs = omodel->mtype;
          printf(" EGADS Info: Ext read failure in %s  %d (EG_loadModel)!\n",
                 name, i);
          break;
        }
        if (ibody > j) {
          omodel->mtype = j;
          mshape->nobjs = omodel->mtype;
          printf(" EGADS Info: Ext read failure in %s  %d %d (EG_loadModel)!\n",
                 name, ibody, j);
          break;
        }
        if (oclass == TESSELLATION) {
          stat = EG_readTess(fp,  mshape->bodies[ibody-1], &mshape->bodies[j]);
          EG_referenceObject(mshape->bodies[ibody-1], mshape->bodies[j]);
        } else {
          stat = EG_readEBody(fp, mshape->bodies[ibody-1], &mshape->bodies[j]);
        }
        if (stat != EGADS_SUCCESS) {
          omodel->mtype = j;
          mshape->nobjs = omodel->mtype;
          printf(" EGADS Info: Ext read failure in %s  %d %d %d (EG_loadModel)!\n",
                 name, oclass, ibody, stat);
          break;
        }
        EG_referenceObject(mshape->bodies[j], omodel);
        EG_removeCntxtRef(mshape->bodies[j]);
        mshape->bodies[j]->topObj = omodel;
      }
    }
  } else {
    printf(" EGADS Info: EGADS Header not found in %s (EG_loadModel)!\n",
           name);
  }

  fclose(fp);
  return EGADS_SUCCESS;
}


void
EG_writeAttr(egAttrs *attrs, FILE *fp)
{
  char c;
  int  namln;

  int    nattr = attrs->nattrs;
  egAttr *attr = attrs->attrs;
  for (int i = 0; i < nattr; i++) {
    if (attr[i].type == ATTRPTR) continue;
    namln = 0;
    if (attr[i].name != NULL) namln = strlen(attr[i].name);
    fprintf(fp, "%d %d %d\n", attr[i].type, namln, attr[i].length);
    if (namln != 0) {
      for (int j = 0; j < namln; j++) {
        c = attr[i].name[j];
        if (c == 32) c = 127;
        fprintf(fp, "%c", c);
      }
      fprintf(fp, "\n");
    }
    if (attr[i].type == ATTRINT) {
      if (attr[i].length == 1) {
        fprintf(fp, "%d\n", attr[i].vals.integer);
      } else {
        for (int j = 0; j < attr[i].length; j++)
          fprintf(fp, "%d ", attr[i].vals.integers[j]);
        fprintf(fp, "\n");
      }
    } else if ((attr[i].type == ATTRREAL) || (attr[i].type == ATTRCSYS)) {
      if (attr[i].length == 1) {
        fprintf(fp, "%19.12le\n", attr[i].vals.real);
      } else {
        for (int j = 0; j < attr[i].length; j++)
          fprintf(fp, "%19.12le ", attr[i].vals.reals[j]);
        fprintf(fp, "\n");
      }
    } else if (attr[i].type == ATTRSTRING) {
      if (attr[i].length != 0) {
        fprintf(fp, "#%s\n", attr[i].vals.string);
      } else {
        /* we should never get here! */
        fprintf(fp, "#\n");
      }
    }
  }
}


int
EG_writeNumAttr(egAttrs *attrs)
{
  int num = 0;
  for (int i = 0; i < attrs->nattrs; i++) {
    if (attrs->attrs[i].type == ATTRPTR) continue;
    num++;
  }

  return num;
}


static void
EG_writeAttrs(const egObject *obj, FILE *fp)
{
  int     i, nsolid, nshell, nface, nloop, nedge, nnode, nattr = 0;
  egAttrs *attrs;

  attrs = (egAttrs *) obj->attrs;
  if (attrs != NULL) nattr = attrs->nattrs;

  if (obj->oclass == MODEL) {

    fprintf(fp, "%d\n", nattr);
    if (nattr != 0) EG_writeAttr(attrs, fp);

  } else {

    egadsBody *pbody = (egadsBody *) obj->blind;
    nnode  = pbody->nodes.map.Extent();
    nedge  = pbody->edges.map.Extent();
    nloop  = pbody->loops.map.Extent();
    nface  = pbody->faces.map.Extent();
    nshell = pbody->shells.map.Extent();
    nsolid = 0;
    if (obj->mtype == SOLIDBODY) nsolid = 1;
    fprintf(fp, "  %d  %d  %d  %d  %d  %d  %d\n", nsolid, nshell, nface,
            nloop, nedge, nnode, nattr);
    if (nattr != 0) EG_writeAttr(attrs, fp);

    for (i = 0; i < nshell; i++) {
      egObject *aobj = pbody->shells.objs[i];
      if (aobj->attrs == NULL) continue;
      attrs = (egAttrs *) aobj->attrs;
      nattr = EG_writeNumAttr(attrs);
      if (nattr <= 0) continue;
      fprintf(fp, "    1 %d %d\n", i, nattr);
      EG_writeAttr(attrs, fp);
    }

    for (i = 0; i < nface; i++) {
      egObject *aobj = pbody->faces.objs[i];
      if (aobj->attrs == NULL) continue;
      attrs = (egAttrs *) aobj->attrs;
      nattr = EG_writeNumAttr(attrs);
      if (nattr <= 0) continue;
      fprintf(fp, "    2 %d %d\n", i, nattr);
      EG_writeAttr(attrs, fp);
    }

    for (i = 0; i < nloop; i++) {
      egObject *aobj = pbody->loops.objs[i];
      if (aobj->attrs == NULL) continue;
      attrs = (egAttrs *) aobj->attrs;
      nattr = EG_writeNumAttr(attrs);
      if (nattr <= 0) continue;
      fprintf(fp, "    3 %d %d\n", i, nattr);
      EG_writeAttr(attrs, fp);
    }

    for (i = 0; i < nedge; i++) {
      egObject *aobj = pbody->edges.objs[i];
      if (aobj->attrs == NULL) continue;
      attrs = (egAttrs *) aobj->attrs;
      nattr = EG_writeNumAttr(attrs);
      if (nattr <= 0) continue;
      fprintf(fp, "    4 %d %d\n", i, nattr);
      EG_writeAttr(attrs, fp);
    }

    for (int i = 0; i < nnode; i++) {
      egObject *aobj = pbody->nodes.objs[i];
      if (aobj->attrs == NULL) continue;
      attrs = (egAttrs *) aobj->attrs;
      nattr = EG_writeNumAttr(attrs);
      if (nattr <= 0) continue;
      fprintf(fp, "    5 %d %d\n", i, nattr);
      EG_writeAttr(attrs, fp);
    }
    fprintf(fp, "    0 0 0\n");

  }

}


static int
EG_writeTess(const egObject *tess, FILE *fp)
{
  int          j, status, len;
  int          ntri, iedge, iface, nnode, nedge, nface, nattr = 0;
  const double *pxyz  = NULL, *puv    = NULL, *pt    = NULL;
  const int    *ptype = NULL, *pindex = NULL, *ptris = NULL, *ptric = NULL;
  egAttrs      *attrs;

  if (tess == NULL)                 return EGADS_NULLOBJ;
  if (tess->magicnumber != MAGIC)   return EGADS_NOTOBJ;
  if (tess->oclass != TESSELLATION) return EGADS_NOTTESS;
  if (tess->blind == NULL)          return EGADS_NODATA;

  // get the body from tessellation
  egTessel  *btess = (egTessel *) tess->blind;
  egObject  *body  = btess->src;

  // get the sizes
  status = EG_getBodyTopos(body, NULL, NODE, &nnode, NULL);
  if (status != EGADS_SUCCESS) return status;
  if (body->oclass == EBODY) {
    status = EG_getBodyTopos(body, NULL, EEDGE, &nedge, NULL);
    if (status != EGADS_SUCCESS) return status;
    status = EG_getBodyTopos(body, NULL, EFACE, &nface, NULL);
    if (status != EGADS_SUCCESS) return status;
  } else {
    status = EG_getBodyTopos(body, NULL, EDGE, &nedge, NULL);
    if (status != EGADS_SUCCESS) return status;
    status = EG_getBodyTopos(body, NULL, FACE, &nface, NULL);
    if (status != EGADS_SUCCESS) return status;
  }

  // write the number of nodes, edges, and faces
  fprintf(fp, " %d %d %d\n", nnode, nedge, nface);

  // write out the edge tessellation
  for (iedge = 0; iedge < nedge; iedge++) {
    status = EG_getTessEdge(tess, iedge+1, &len, &pxyz, &pt);
    if (status != EGADS_SUCCESS) return status;
    fprintf(fp, " %d\n", len);
    if (len == 0) continue;
    for (j = 0; j < len; j++)
      fprintf(fp, "%20.16le %20.16le %20.16le %20.16le\n", pxyz[3*j],
              pxyz[3*j+1], pxyz[3*j+2], pt[j]);
  }

  // write out face tessellations
  for (iface = 0; iface < nface; iface++) {
    status = EG_getTessFace(tess, iface+1, &len, &pxyz, &puv, &ptype, &pindex,
                            &ntri, &ptris, &ptric);
    if ((status != EGADS_SUCCESS) && (status != EGADS_NODATA)) return status;
    fprintf(fp, " %d %d\n", len, ntri);
    if ((len == 0) || (ntri == 0)) continue;
    for (j = 0; j < len; j++)
      fprintf(fp, "%20.16le %20.16le %20.16le %20.16le %20.16le %d %d\n",
              pxyz[3*j], pxyz[3*j+1], pxyz[3*j+2], puv[2*j], puv[2*j+1],
              ptype[j], pindex[j]);
    for (j = 0; j < ntri; j++)
      fprintf(fp, "%d %d %d %d %d %d\n", ptris[3*j], ptris[3*j+1], ptris[3*j+2],
              ptric[3*j], ptric[3*j+1], ptric[3*j+2]);
  }

  // write out the tessellation attributes
  attrs = (egAttrs *) tess->attrs;
  if (attrs != NULL) nattr = attrs->nattrs;
  fprintf(fp, "%d\n", nattr);
  if (nattr != 0) EG_writeAttr(attrs, fp);

  return EGADS_SUCCESS;
}


//Based on STEPCAFControl_Writer FindEntities
static Standard_Integer
EG_FindEntities(const Handle(Transfer_FinderProcess)& theFP,
                const TopoDS_Shape& theShape,
                TopLoc_Location& theLocation,
                TColStd_SequenceOfTransient& theSeqRI)
{
  Handle(StepRepr_RepresentationItem)
    anItem = STEPConstruct::FindEntity(theFP, theShape, theLocation);

  if (!anItem.IsNull())
  {
    theSeqRI.Append(anItem);
    return 1;
  }

  // may be S was splited during shape processing
  Handle(TransferBRep_ShapeMapper)
    aMapper = TransferBRep::ShapeMapper(theFP, theShape);
  Handle(Transfer_Binder) aBinder = theFP->Find(aMapper);
  if (aBinder.IsNull()) return 0;

  Handle(Transfer_TransientListBinder) aTransientListBinder =
    Handle(Transfer_TransientListBinder)::DownCast(aBinder);
  Standard_Integer aResCount = 0;
  if (aTransientListBinder.IsNull() &&
      theShape.ShapeType() == TopAbs_COMPOUND) {
    for (TopoDS_Iterator anIter(theShape); anIter.More(); anIter.Next()) {
      Handle(StepRepr_RepresentationItem) aLocalItem =
        STEPConstruct::FindEntity(theFP, anIter.Value(), theLocation);
      if (aLocalItem.IsNull()) continue;
      aResCount++;
      theSeqRI.Append(aLocalItem);
    }
  } else if (!aTransientListBinder.IsNull()) {
    const Standard_Integer aNbTransient = aTransientListBinder->NbTransients();
    for (Standard_Integer anInd = 1; anInd <= aNbTransient; anInd++) {
      Handle(Standard_Transient)
        anEntity = aTransientListBinder->Transient(anInd);
      anItem = Handle(StepRepr_RepresentationItem)::DownCast(anEntity);
      if (anItem.IsNull()) continue;
      aResCount++;
      theSeqRI.Append(anItem);
    }
  }
  return aResCount;
}


//Based on STEPConstruct_Styles::EncodeColor
Handle(StepVisual_Colour)
EG_STEPEncodeColor(const Quantity_Color &C,
                   STEPConstruct_DataMapOfAsciiStringTransient &DPDCs,
                   STEPConstruct_DataMapOfPointTransient &ColRGBs)
{
  // detect if color corresponds to one of pre-defined colors
  Standard_CString cName = 0;
       if ( C == Quantity_Color(Quantity_NOC_RED)                ) cName = "red";
  else if ( C == Quantity_Color(Quantity_NOC_GREEN)              ) cName = "green";
  else if ( C == Quantity_Color(Quantity_NOC_BLUE1)              ) cName = "blue";
  else if ( C == Quantity_Color(Quantity_NOC_YELLOW)             ) cName = "yellow";
  else if ( C == Quantity_Color(Quantity_NOC_MAGENTA1)           ) cName = "magenta";
  else if ( C == Quantity_Color(Quantity_NOC_CYAN1)              ) cName = "cyan";
  else if ( C == Quantity_Color(Quantity_NOC_BLACK)              ) cName = "black";
  else if ( C == Quantity_Color(Quantity_NOC_WHITE)              ) cName = "white";
  else if ( C == Quantity_Color(1.0, 0.5, 0.5, Quantity_TOC_RGB) ) cName = "lred";
  else if ( C == Quantity_Color(0.5, 1.0, 0.5, Quantity_TOC_RGB) ) cName = "lgreen";
  else if ( C == Quantity_Color(0.5, 0.5, 1.0, Quantity_TOC_RGB) ) cName = "lblue";

  if (cName) {
    Handle(StepVisual_DraughtingPreDefinedColour) ColPr;
    TCollection_AsciiString aName(cName);
    if(DPDCs.IsBound(aName)) {
      ColPr = Handle(StepVisual_DraughtingPreDefinedColour)::DownCast(DPDCs.Find(aName));
      if(!ColPr.IsNull()) return ColPr;
    }
    ColPr = new StepVisual_DraughtingPreDefinedColour;
    Handle(StepVisual_PreDefinedItem) preDef = new StepVisual_PreDefinedItem;
    preDef->Init(new TCollection_HAsciiString(cName));
    ColPr->SetPreDefinedItem(preDef);
    DPDCs.Bind(aName,ColPr);
    return ColPr;
  } else {
    Handle(StepVisual_ColourRgb) ColRGB;
    gp_Pnt P;
    C.Values (P.ChangeCoord().ChangeData()[0],
              P.ChangeCoord().ChangeData()[1],
              P.ChangeCoord().ChangeData()[2],
              Quantity_TOC_RGB);
    if(ColRGBs.IsBound(P)) {
      ColRGB = Handle(StepVisual_ColourRgb)::DownCast(ColRGBs.Find(P));
      if(!ColRGB.IsNull()) return ColRGB;
    }
    Handle(TCollection_HAsciiString) ColName = new TCollection_HAsciiString("");
    ColRGB = new StepVisual_ColourRgb;
    ColRGB->Init(ColName, P.Coord(1), P.Coord(2), P.Coord(3));
    ColRGBs.Bind(P,ColRGB);
    return ColRGB;
  }
}


//Based on STEPCAFControl_Writer::writeColors
//     and STEPCAFControl_Writer MakeSTEPStyles
static void
EG_stepSetColors(const Quantity_Color* aSurfCol,
                 const Quantity_Color* aCurveCol,
                 const Handle(XSControl_WorkSession)& theWS,
                 MoniTool_DataMapOfShapeTransient& theMapCompMDGPR,
                 const TopoDS_Shape& aTopSh,
                 const TopoDS_Shape& aShape)
{

  STEPConstruct_Styles Styles(theWS);
  STEPConstruct_DataMapOfAsciiStringTransient DPDCs;
  STEPConstruct_DataMapOfPointTransient ColRGBs;

  Handle(StepRepr_RepresentationContext) aContext = Styles.FindContext(aTopSh);
  if (aContext.IsNull()) return;

  // iterate on subshapes and create STEP styles
  Handle(StepVisual_StyledItem) anOverride;
  TopTools_MapOfShape aMap;

  // translate colors to STEP
  Handle(StepVisual_Colour) aSurfColor, aCurvColor;
  if (aSurfCol != NULL)
    aSurfColor = EG_STEPEncodeColor(*aSurfCol, DPDCs, ColRGBs);
  if (aCurveCol != NULL)
    aCurvColor = EG_STEPEncodeColor(*aCurveCol, DPDCs, ColRGBs);


  TopLoc_Location aLocation;
  TColStd_SequenceOfTransient aSeqRI;
  Standard_Integer aNbEntities = EG_FindEntities(Styles.FinderProcess(),
                                                 aShape, aLocation, aSeqRI);
  if (aNbEntities <= 0)
    std::cout << "EGADS Warning: Cannot find RI for " <<
                 aShape.TShape()->DynamicType()->Name() << std::endl;

  for (TColStd_SequenceOfTransient::Iterator anEntIter(aSeqRI);
      anEntIter.More(); anEntIter.Next()) {
    const Handle(StepRepr_RepresentationItem)& anItem =
        Handle(StepRepr_RepresentationItem)::DownCast(anEntIter.Value());
    Handle(StepVisual_PresentationStyleAssignment) aPSA;
#if CASVER < 760
    aPSA = Styles.MakeColorPSA(anItem, aSurfColor, aCurvColor, Standard_False);
#else
    Standard_Real aRenderTransp = 0.0;
    aPSA = Styles.MakeColorPSA(anItem, aSurfColor, aCurvColor, aSurfColor,
                               aRenderTransp, Standard_False);
#endif
    Handle(StepVisual_StyledItem) theOverride;
    Styles.AddStyle(anItem, aPSA, theOverride);
  }

  const Handle(XSControl_TransferWriter)& aTW = theWS->TransferWriter();
  const Handle(Transfer_FinderProcess)& aFP = aTW->FinderProcess();
  Handle(StepData_StepModel)
    aStepModel = Handle(StepData_StepModel)::DownCast(aFP->Model());

  // create MDGPR and record it in model
  Handle(StepVisual_MechanicalDesignGeometricPresentationRepresentation) aMDGPR;

  if (!theMapCompMDGPR.IsBound(aTopSh)) {
#if CASVER < 780
    Styles.CreateMDGPR(aContext, aMDGPR);
#else
    Styles.CreateMDGPR(aContext, aMDGPR, aStepModel);
#endif
    if (!aMDGPR.IsNull()) theMapCompMDGPR.Bind(aTopSh, aMDGPR);
  } else {
    aMDGPR = Handle(StepVisual_MechanicalDesignGeometricPresentationRepresentation)::DownCast(theMapCompMDGPR.Find(aTopSh));
  }

  if (aMDGPR.IsNull()) return;

  Handle(StepRepr_HArray1OfRepresentationItem) anOldItems = aMDGPR->Items();
  Standard_Integer oldLengthlen = 0;
  if (!anOldItems.IsNull())
    oldLengthlen = anOldItems->Length();
  const Standard_Integer aNbIt = oldLengthlen + Styles.NbStyles();
  if (!aNbIt) return;
  Handle(StepRepr_HArray1OfRepresentationItem) aNewItems =
      new StepRepr_HArray1OfRepresentationItem(1, aNbIt);
  Standard_Integer anElemInd = 1;
  for (Standard_Integer aStyleInd = 1; aStyleInd <= oldLengthlen; aStyleInd++) {
    aNewItems->SetValue(anElemInd++, anOldItems->Value(aStyleInd));
  }
  for (Standard_Integer aStyleInd = 1; aStyleInd <= Styles.NbStyles(); aStyleInd++) {
    aNewItems->SetValue(anElemInd++, Styles.Style(aStyleInd));
  }

  if (aNewItems->Length() > 0) aMDGPR->SetItems(aNewItems);
}


static bool
EG_quantityColor(const int aType,
                 const int aLen,
                 const char *str,
                 const double *reals,
                 Quantity_Color& aCol)
{
  bool found = false;
  if ((aType == ATTRREAL && aLen == 3)) {
    aCol = Quantity_Color(reals[0], reals[1], reals[2], Quantity_TOC_RGB);
    found = true;
  } else if (aType == ATTRSTRING) {
    if        (strcasecmp(str, "red"    ) == 0) {
      aCol = Quantity_Color(Quantity_NOC_RED1    ); found = true;
    } else if (strcasecmp(str, "green"  ) == 0) {
      aCol = Quantity_Color(Quantity_NOC_GREEN1  ); found = true;
    } else if (strcasecmp(str, "blue"   ) == 0) {
      aCol = Quantity_Color(Quantity_NOC_BLUE1   ); found = true;
    } else if (strcasecmp(str, "yellow" ) == 0) {
      aCol = Quantity_Color(Quantity_NOC_YELLOW1 ); found = true;
    } else if (strcasecmp(str, "magenta") == 0) {
      aCol = Quantity_Color(Quantity_NOC_MAGENTA1); found = true;
    } else if (strcasecmp(str, "cyan"   ) == 0) {
      aCol = Quantity_Color(Quantity_NOC_CYAN1   ); found = true;
    } else if (strcasecmp(str, "white"  ) == 0) {
      aCol = Quantity_Color(Quantity_NOC_WHITE   ); found = true;
    } else if (strcasecmp(str, "black"  ) == 0) {
      aCol = Quantity_Color(Quantity_NOC_BLACK   ); found = true;
    } else if (strcasecmp(str, "lred"   ) == 0) {
      aCol = Quantity_Color(1.0, 0.5, 0.5, Quantity_TOC_RGB); found = true;
    } else if (strcasecmp(str, "lgreen" ) == 0) {
      aCol = Quantity_Color(0.5, 1.0, 0.5, Quantity_TOC_RGB); found = true;
    } else if (strcasecmp(str, "lblue"  ) == 0) {
      aCol = Quantity_Color(0.5, 0.5, 1.0, Quantity_TOC_RGB); found = true;
    }
  }

  return found;
}

static Handle(IGESGraph_Color)
EG_igesGraphColor(const Quantity_Color& color)
{
  Handle(IGESGraph_Color) colorEntity = new IGESGraph_Color;

  const double rgb[3] = {color.Red(), color.Green(), color.Blue()};

         if (rgb[0] == 1   && rgb[1] == 0   && rgb[2] == 0  ) {
    colorEntity->Init(100,   0,   0, new TCollection_HAsciiString("red"    ));
  } else if (rgb[0] == 0   && rgb[1] == 1   && rgb[2] == 0  ) {
    colorEntity->Init(  0, 100,   0, new TCollection_HAsciiString("green"  ));
  } else if (rgb[0] == 0   && rgb[1] == 0   && rgb[2] == 1  ) {
    colorEntity->Init(  0,   0, 100, new TCollection_HAsciiString("blue"   ));
  } else if (rgb[0] == 1   && rgb[1] == 1   && rgb[2] == 0  ) {
    colorEntity->Init(100, 100,   0, new TCollection_HAsciiString("yellow" ));
  } else if (rgb[0] == 1   && rgb[1] == 0   && rgb[2] == 1  ) {
    colorEntity->Init(100,   0, 100, new TCollection_HAsciiString("magenta"));
  } else if (rgb[0] == 0   && rgb[1] == 1   && rgb[2] == 1  ) {
    colorEntity->Init(  0, 100, 100, new TCollection_HAsciiString("cyan"   ));
  } else if (rgb[0] == 1   && rgb[1] == 1   && rgb[2] == 1  ) {
    colorEntity->Init(100, 100, 100, new TCollection_HAsciiString("white"  ));
  } else if (rgb[0] == 0   && rgb[1] == 0   && rgb[2] == 0  ) {
    colorEntity->Init(  0,   0,   0, new TCollection_HAsciiString("black"  ));
  } else if (rgb[0] == 1   && rgb[1] == 0.5 && rgb[2] == 0.5) {
    colorEntity->Init(100,  50,  50, new TCollection_HAsciiString("lred"   ));
  } else if (rgb[0] == 0.5 && rgb[1] == 1   && rgb[2] == 0.5) {
    colorEntity->Init( 50, 100,  50, new TCollection_HAsciiString("lgreen" ));
  } else if (rgb[0] == 0.5 && rgb[1] == 0.5 && rgb[2] == 1  ) {
    colorEntity->Init( 50,  50, 100, new TCollection_HAsciiString("lblue"  ));
  } else {
    colorEntity->Init(100*rgb[0], 100*rgb[1], 100*rgb[2], NULL);
  }

  return colorEntity;
}


static void
EG_setSTEPprops(const Handle(XSControl_WorkSession) &WS, int nbody,
                const egObject **bodies, Handle(Transfer_FinderProcess) FP)
{
  int          i, j, stat, aType, aLen, nshell, nface, nloop, nedge, nnode;
  const int    *ints;
  const double *reals;
  const char   *str;
  TopLoc_Location Loc;
  Handle(StepRepr_RepresentationItem) r;
  MoniTool_DataMapOfShapeTransient theMapCompMDGPR;
  STEPConstruct_Styles Styles(WS);
  const Handle(Interface_InterfaceModel) &aModel = WS->Model();

#ifdef WRITECSYS
  // Get working data & place to store Body CSYSs
  Handle(StepRepr_HArray1OfRepresentationItem) reprItems =
    new StepRepr_HArray1OfRepresentationItem;
#endif

  for (j = 0; j < nbody; j++) {
    egadsBody *pbody = (egadsBody *) bodies[j]->blind;
    nnode  = pbody->nodes.map.Extent();
    nedge  = pbody->edges.map.Extent();
    nloop  = pbody->loops.map.Extent();
    nface  = pbody->faces.map.Extent();
    nshell = pbody->shells.map.Extent();

    r = STEPConstruct::FindEntity(FP, pbody->shape, Loc);
    if (bodies[j]->mtype == WIREBODY) {
      Handle(StepShape_GeometricCurveSet)
        curveSet = Handle(StepShape_GeometricCurveSet)::DownCast(r);
      if (!curveSet.IsNull()) {
        for (i = 0; i < curveSet->NbElements(); i++) {
          egObject  *aobj  = pbody->edges.objs[i];
          if (aobj  == NULL) continue;
          egadsEdge *pedge = (egadsEdge *) aobj->blind;
          if (pedge == NULL) continue;
          stat = EG_attributeRet(aobj, "Name", &aType, &aLen, &ints, &reals, &str);
          if (stat == EGADS_SUCCESS && aType == ATTRSTRING) {
            StepShape_GeometricSetSelect curveSelect = curveSet->ElementsValue(i+1);
            Handle(StepGeom_Curve) curve = curveSelect.Curve();
            if (curve.IsNull()) continue;
            Handle(TCollection_HAsciiString) id = new TCollection_HAsciiString(str);
            curve->SetName(id);

            // If it's a trimmed curve, set name on basis curve as well
            Handle(StepGeom_TrimmedCurve)
              trimcurve = Handle(StepGeom_TrimmedCurve)::DownCast(curve);
            if (!trimcurve.IsNull()) {
              Handle(StepGeom_Curve) basiscurve = trimcurve->BasisCurve();
              if (!basiscurve.IsNull()) {
                basiscurve->SetName(id);
              }
            }
          }

#ifdef STEP_WIRE_BODY_EDGE_COLORS_DOES_NOT_WORK
          stat = EG_attributeRet(aobj, "Color", &aType, &aLen, &ints, &reals, &str);
          if (stat == EGADS_SUCCESS) {
            Quantity_Color aCurvCol;
            if (EG_quantityColor(aType, aLen, str, reals, aCurvCol))
              EG_stepSetColors(NULL, &aCurvCol,
                               WS, theMapCompMDGPR,
                               pbody->shape, r);
          }
#endif
        }
      }
    }

    stat   = EG_attributeRet(bodies[j], "Name", &aType, &aLen, &ints, &reals,
                             &str);
    if (stat == EGADS_SUCCESS && aType == ATTRSTRING) {
      if (!r.IsNull()) {
        Handle(TCollection_HAsciiString) id = new TCollection_HAsciiString(str);
        r->SetName(id);
      }
    }

    bool ColorEdge = false, ColorFace = false;
    Quantity_Color aCurvCol, aSurfCol;
    stat = EG_attributeRet(bodies[j], "ColorEdge", &aType, &aLen, &ints, &reals,
                           &str);
    if (stat == EGADS_SUCCESS) {
      ColorEdge = EG_quantityColor(aType, aLen, str, reals, aCurvCol);
    }
    stat = EG_attributeRet(bodies[j], "ColorFace", &aType, &aLen, &ints, &reals,
                           &str);
    if (stat == EGADS_SUCCESS) {
      ColorFace = EG_quantityColor(aType, aLen, str, reals, aSurfCol);
    }
    if (ColorEdge || ColorFace)
      EG_stepSetColors(ColorFace ? &aSurfCol : NULL, ColorEdge ? &aCurvCol : NULL,
                       WS, theMapCompMDGPR,
                       pbody->shape, pbody->shape);

#ifdef WRITECSYS
    // find unattached Body CSYSs
    int nAttr;
    stat = EG_attributeNum(bodies[j], &nAttr);
    if (stat == EGADS_SUCCESS)
      for (i = 1; i <= nAttr; i++) {
        const char *name;
        stat = EG_attributeGet(bodies[j], i, &name, &aType, &aLen, &ints,
                               &reals, &str);
        if (stat == EGADS_SUCCESS)
          if ((aType == ATTRCSYS) && (aLen == 9)) {
            gp_Pnt pntc(reals[ 9], reals[10], reals[11]);
            gp_Dir dirx(reals[12], reals[13], reals[14]);
            gp_Dir diry(reals[15], reals[16], reals[17]);
            gp_Ax2 aDTAxis(pntc, dirx, diry);
            GeomToStep_MakeAxis2Placement3d anAxisMaker(aDTAxis);
            Handle(StepGeom_Axis2Placement3d) anA2P3D = anAxisMaker.Value();
            anA2P3D->SetName(new TCollection_HAsciiString(name));
            int len = reprItems->Length() + 1;
            reprItems->Resize(1, len, Standard_True);
            reprItems->SetValue(len, anA2P3D);
          }
      }
#endif

    for (i = 0; i < nshell; i++) {
      egObject   *aobj   = pbody->shells.objs[i];
      if (aobj   == NULL) continue;
      egadsShell *pshell = (egadsShell *) aobj->blind;
      if (pshell == NULL) continue;
      stat = EG_attributeRet(aobj, "Name", &aType, &aLen, &ints, &reals, &str);
      if (stat  != EGADS_SUCCESS) continue;
      if (aType != ATTRSTRING)    continue;
      r = STEPConstruct::FindEntity(FP, pshell->shell, Loc);
      if (r.IsNull()) continue;
      Handle(StepShape_ShellBasedSurfaceModel)
        shellModel = Handle(StepShape_ShellBasedSurfaceModel)::DownCast(r);
      if (shellModel.IsNull()) continue;
      Handle(TCollection_HAsciiString) id = new TCollection_HAsciiString(str);
      r->SetName(id);
      for (Standard_Integer i = 1; i <= shellModel->NbSbsmBoundary(); i++) {
        StepShape_Shell shell = shellModel->SbsmBoundaryValue(i);
        if (!shell.OpenShell().IsNull()  ) shell.OpenShell()->SetName(id);
        if (!shell.ClosedShell().IsNull()) shell.ClosedShell()->SetName(id);
      }
    }

    for (i = 0; i < nface; i++) {
      egObject  *aobj  = pbody->faces.objs[i];
      if (aobj  == NULL) continue;
      egadsFace *pface = (egadsFace *) aobj->blind;
      if (pface == NULL) continue;
      stat = EG_attributeRet(aobj, "Name", &aType, &aLen, &ints, &reals, &str);
      if (stat == EGADS_SUCCESS && aType == ATTRSTRING) {
        r = STEPConstruct::FindEntity(FP, pface->face, Loc);
        if (!r.IsNull()) {
          Handle(TCollection_HAsciiString) id = new TCollection_HAsciiString(str);
          r->SetName(id);

          // Set the name on the underlying surface as well
          Handle(StepShape_FaceSurface)
            faceSurf = Handle(StepShape_FaceSurface)::DownCast(r);
          if (!faceSurf.IsNull()) {
            Handle(StepGeom_Surface) surface = faceSurf->FaceGeometry();
            if (!surface.IsNull()) {
              surface->SetName(id);
            }
          }
        }
      }

      stat = EG_attributeRet(aobj, "Color", &aType, &aLen, &ints, &reals, &str);
      if (stat == EGADS_SUCCESS) {
        Quantity_Color aSurfCol;
        if (EG_quantityColor(aType, aLen, str, reals, aSurfCol))
          EG_stepSetColors(&aSurfCol, NULL,
                           WS, theMapCompMDGPR,
                           pbody->shape, pface->face);
      }
    }

    for (i = 0; i < nloop; i++) {
      egObject *aobj   = pbody->loops.objs[i];
      if (aobj  == NULL) continue;
      egadsLoop *ploop = (egadsLoop *) aobj->blind;
      if (ploop == NULL) continue;
      stat = EG_attributeRet(aobj, "Name", &aType, &aLen, &ints, &reals, &str);
      if (stat  != EGADS_SUCCESS) continue;
      if (aType != ATTRSTRING)    continue;
      r = STEPConstruct::FindEntity(FP, ploop->loop, Loc);
      if (r.IsNull()) continue;
      Handle(StepShape_GeometricCurveSet)
        curveSet = Handle(StepShape_GeometricCurveSet)::DownCast(r);
      if (curveSet.IsNull()) continue;
      Handle(TCollection_HAsciiString) id = new TCollection_HAsciiString(str);
      r->SetName(id);
      for (Standard_Integer i = 1; i <= curveSet->NbElements(); i++) {
        StepShape_GeometricSetSelect curveSelect = curveSet->ElementsValue(i);
        Handle(StepGeom_Curve) curve = curveSelect.Curve();
        if (curve.IsNull()) continue;
        curve->SetName(id);
        
        // If it's a trimmed curve, set name on basis curve as well
        Handle(StepGeom_TrimmedCurve)
          trimcurve = Handle(StepGeom_TrimmedCurve)::DownCast(curve);
        if (!trimcurve.IsNull()) {
          Handle(StepGeom_Curve) basiscurve = trimcurve->BasisCurve();
          if (!basiscurve.IsNull()) {
            basiscurve->SetName(id);
          }
        }

        // If it's a surface curve, set name on 3d curve as well
        Handle(StepGeom_SurfaceCurve)
          surfcurve = Handle(StepGeom_SurfaceCurve)::DownCast(curve);
        if (!surfcurve.IsNull()) {
          Handle(StepGeom_Curve) curve3d = surfcurve->Curve3d();
          if (!curve3d.IsNull()) {
            curve3d->SetName(id);
          }
        }
      }
    }

    for (i = 0; i < nedge; i++) {
      egObject  *aobj  = pbody->edges.objs[i];
      if (aobj  == NULL) continue;
      egadsEdge *pedge = (egadsEdge *) aobj->blind;
      if (pedge == NULL) continue;
      stat = EG_attributeRet(aobj, "Name", &aType, &aLen, &ints, &reals, &str);
      if (stat == EGADS_SUCCESS && aType == ATTRSTRING) {
        r = STEPConstruct::FindEntity(FP, pedge->edge, Loc);
        if (r.IsNull()) continue;
        Handle(TCollection_HAsciiString) id = new TCollection_HAsciiString(str);
        r->SetName(id);

        // Set the name on the underlying curve as well
        Handle(StepShape_EdgeCurve)
          edgeCurve = Handle(StepShape_EdgeCurve)::DownCast(r);
        if (!edgeCurve.IsNull()) {
          Handle(StepGeom_Curve) curve = edgeCurve->EdgeGeometry();
          if (!curve.IsNull()) {
            curve->SetName(id);

            // If it's a trimmed curve, set name on basis curve as well
            Handle(StepGeom_TrimmedCurve)
              trimcurve = Handle(StepGeom_TrimmedCurve)::DownCast(curve);
            if (!trimcurve.IsNull()) {
              Handle(StepGeom_Curve) basiscurve = trimcurve->BasisCurve();
              if (!basiscurve.IsNull()) {
                basiscurve->SetName(id);
              }
            }

            // If it's a surface curve, set name on 3d curve as well
            Handle(StepGeom_SurfaceCurve)
              surfcurve = Handle(StepGeom_SurfaceCurve)::DownCast(curve);
            if (!surfcurve.IsNull()) {
              Handle(StepGeom_Curve) curve3d = surfcurve->Curve3d();
              if (!curve3d.IsNull()) {
                curve3d->SetName(id);
              }
            }
          }
        }
      }

      stat = EG_attributeRet(aobj, "Color", &aType, &aLen, &ints, &reals, &str);
      if (stat == EGADS_SUCCESS) {
        Quantity_Color aCurveCol;
        if (EG_quantityColor(aType, aLen, str, reals, aCurveCol))
          EG_stepSetColors(NULL, &aCurveCol,
                           WS, theMapCompMDGPR,
                           pbody->shape, pedge->edge);
      }
    }

    for (int i = 0; i < nnode; i++) {
      egObject  *aobj  = pbody->nodes.objs[i];
      if (aobj  == NULL) continue;
      egadsNode *pnode = (egadsNode *) aobj->blind;
      if (pnode == NULL) continue;
      stat = EG_attributeRet(aobj, "Name", &aType, &aLen, &ints, &reals, &str);
      if (stat  != EGADS_SUCCESS) continue;
      if (aType != ATTRSTRING)    continue;
      r = STEPConstruct::FindEntity(FP, pnode->node, Loc);
      if (r.IsNull()) continue;
      Handle(TCollection_HAsciiString) id = new TCollection_HAsciiString(str);
      r->SetName(id);

      // Set the name on the underlying point as well
      Handle(StepShape_VertexPoint)
        vertPoint = Handle(StepShape_VertexPoint)::DownCast(r);
      if (!vertPoint.IsNull()) {
        Handle(StepGeom_Point) point = vertPoint->VertexGeometry();
        if (!point.IsNull()) {
          point->SetName(id);
        }
      }
    }
  }

#ifdef WRITECSYS
  // do we need to add any Body CSYSs?
  if (reprItems->Length() > 0) {
    Handle(StepRepr_RepresentationContext) dummyRC;
    Handle(StepRepr_ConstructiveGeometryRepresentation) theRepr =
      new StepRepr_ConstructiveGeometryRepresentation();
    theRepr->Init(new TCollection_HAsciiString("supplemental geometry"),
                  reprItems, dummyRC);
    aModel->AddWithRefs(theRepr);
  }
#endif

  // register all MDGPRs in model
  for (MoniTool_DataMapIteratorOfDataMapOfShapeTransient anItr(theMapCompMDGPR);
       anItr.More(); anItr.Next())
  {
    aModel->AddWithRefs(anItr.Value());
  }
}


static void
EG_getShapeData(int nbody, const egObject **bodies,
                egadsShapeData &shapeData)
{
  int          i, j, stat, aType, aLen, nshell, nface, nloop, nedge, nnode;
  const int    *ints;
  const double *reals;
  const char   *str;

  for (j = 0; j < nbody; j++) {
    egadsBody *pbody = (egadsBody *) bodies[j]->blind;
    nnode  = pbody->nodes.map.Extent();
    nedge  = pbody->edges.map.Extent();
    nloop  = pbody->loops.map.Extent();
    nface  = pbody->faces.map.Extent();
    nshell = pbody->shells.map.Extent();

    stat   = EG_attributeRet(bodies[j], "Name", &aType, &aLen, &ints, &reals,
                             &str);
    if (stat == EGADS_SUCCESS && aType == ATTRSTRING) {
      shapeData.Add(pbody->shape, str);
    }

#ifdef IGES_DOES_NOT_SUPPORT_BODY_COLORS
    if (bodies[j]->mtype == WIREBODY) {
      stat = EG_attributeRet(bodies[j], "ColorEdge", &aType, &aLen, &ints,
                             &reals, &str);
    } else {
      stat = EG_attributeRet(bodies[j], "ColorFace", &aType, &aLen, &ints,
                             &reals, &str);
    }
    if (stat == EGADS_SUCCESS) {
      Quantity_Color aCol;
      if (EG_quantityColor(aType, aLen, str, reals, aCol)) {
        shapeData.Add(pbody->shape, aCol);
      }
    }
#endif

    for (i = 0; i < nshell; i++) {
      egObject   *aobj   = pbody->shells.objs[i];
      if (aobj   == NULL) continue;
      egadsShell *pshell = (egadsShell *) aobj->blind;
      if (pshell == NULL) continue;
      stat = EG_attributeRet(aobj, "Name", &aType, &aLen, &ints, &reals, &str);
      if (stat  != EGADS_SUCCESS) continue;
      if (aType != ATTRSTRING)    continue;
      shapeData.Add(pshell->shell, str);
    }

    for (i = 0; i < nface; i++) {
      egObject  *aobj  = pbody->faces.objs[i];
      if (aobj  == NULL) continue;
      egadsFace *pface = (egadsFace *) aobj->blind;
      if (pface == NULL) continue;
      stat = EG_attributeRet(aobj, "Name", &aType, &aLen, &ints, &reals, &str);
      if (stat == EGADS_SUCCESS && aType == ATTRSTRING) {
        shapeData.Add(pface->face, str);
      }

      stat = EG_attributeRet(aobj, "Color", &aType, &aLen, &ints, &reals, &str);
      if (stat == EGADS_SUCCESS) {
        Quantity_Color aSurfCol;
        if (EG_quantityColor(aType, aLen, str, reals, aSurfCol) ) {
          shapeData.Add(pface->face, aSurfCol);
        }
      }
    }

    for (i = 0; i < nloop; i++) {
      egObject *aobj   = pbody->loops.objs[i];
      if (aobj  == NULL) continue;
      egadsLoop *ploop = (egadsLoop *) aobj->blind;
      if (ploop == NULL) continue;
      stat = EG_attributeRet(aobj, "Name", &aType, &aLen, &ints, &reals, &str);
      if (stat  != EGADS_SUCCESS) continue;
      if (aType != ATTRSTRING)    continue;
      shapeData.Add(ploop->loop, str);
    }

    for (i = 0; i < nedge; i++) {
      egObject  *aobj  = pbody->edges.objs[i];
      if (aobj  == NULL) continue;
      egadsEdge *pedge = (egadsEdge *) aobj->blind;
      if (pedge == NULL) continue;

      stat = EG_attributeRet(aobj, "Name", &aType, &aLen, &ints, &reals, &str);
      if (stat == EGADS_SUCCESS && aType == ATTRSTRING) {
        shapeData.Add(pedge->edge, str);
      }

      stat = EG_attributeRet(aobj, "Color", &aType, &aLen, &ints, &reals, &str);
      if (stat == EGADS_SUCCESS) {
        Quantity_Color aCurvCol;
        if (EG_quantityColor(aType, aLen, str, reals, aCurvCol) ) {
          shapeData.Add(pedge->edge, aCurvCol);
        }
      }
    }

    for (int i = 0; i < nnode; i++) {
      egObject  *aobj  = pbody->nodes.objs[i];
      if (aobj  == NULL) continue;
      egadsNode *pnode = (egadsNode *) aobj->blind;
      if (pnode == NULL) continue;
      stat = EG_attributeRet(aobj, "Name", &aType, &aLen, &ints, &reals, &str);
      if (stat  != EGADS_SUCCESS) continue;
      if (aType != ATTRSTRING)    continue;
      shapeData.Add(pnode->node, str);
    }
  }
}


static void
EG_setIGESprop(const egadsShapeData &shapeData,
               IGESControl_Writer &iWrite)
{
  int i;
  Handle(IGESData_IGESEntity) ent;
  const Handle(Transfer_FinderProcess) &FP = iWrite.TransferProcess();

#if CASVER >= 780
  NCollection_DataMap<Quantity_Color,Handle(IGESGraph_Color)> colors;
#endif
  NCollection_DataMap<TCollection_AsciiString,Handle(IGESBasic_Name)> names;

  if (shapeData.LabelExtent() > 0) {
    for (i = 1; i <= shapeData.LabelExtent(); i++) {
      const char *str = shapeData.label(i).shapeName;
      if (str == NULL) continue;
      TopoDS_Shape aShape = shapeData.LabelFindKey(i);

      Handle(TransferBRep_ShapeMapper)
        mapper = TransferBRep::ShapeMapper(FP, aShape);
      if (FP->FindTypedTransient(mapper, STANDARD_TYPE(IGESData_IGESEntity), ent)) {
        // Set short name (max 8 chars)
        // r->SetLabel(new TCollection_HAsciiString("SHRTNAME"));

        //Add long name
        Handle(IGESBasic_Name) nameEntity;
        if (names.IsBound(str)) {
          nameEntity = names.Find(str);
        } else {
          nameEntity = new IGESBasic_Name();
          nameEntity->Init(1, new TCollection_HAsciiString(str));
          names.Bind(str, nameEntity);
          iWrite.Model()->AddEntity(nameEntity);
        }
        ent->AddProperty(nameEntity);
      }
    }
  }

  if (shapeData.ColorExtent() > 0) {
    for (i = 1; i <= shapeData.ColorExtent(); i++) {
      const Quantity_Color& color = shapeData.color(i);
      TopoDS_Shape aShape = shapeData.ColorFindKey(i);

      Handle(TransferBRep_ShapeMapper)
        mapper = TransferBRep::ShapeMapper(FP, aShape);
      Handle(TransferBRep_ShapeMapper)
        aNoLocMapper = TransferBRep::ShapeMapper(FP, aShape.Located(TopLoc_Location()));
      if (FP->FindTypedTransient(mapper, STANDARD_TYPE(IGESData_IGESEntity), ent) ||
          FP->FindTypedTransient(aNoLocMapper, STANDARD_TYPE(IGESData_IGESEntity), ent)) {

        Handle(IGESGraph_Color) colorEntity;
#if CASVER >= 780
        if (colors.IsBound(color)) {
          colorEntity = colors.Find(color);
        } else {
          colorEntity =  EG_igesGraphColor(color);
          colors.Bind(color, colorEntity);
          iWrite.Model()->AddWithRefs(colorEntity,
                                      IGESSelect_WorkLibrary::DefineProtocol());
        }
#else
        colorEntity =  EG_igesGraphColor(color);
        iWrite.Model()->AddWithRefs(colorEntity,
                                    IGESSelect_WorkLibrary::DefineProtocol());
#endif

        Standard_Integer rank = EG_igesEncodeColor ( color );

        ent->InitColor ( colorEntity, rank );
        Handle(IGESSolid_Face) ent_f = Handle(IGESSolid_Face)::DownCast(ent);
        if (!ent_f.IsNull())
          if (!ent_f->Surface().IsNull())
            ent_f->Surface()->InitColor ( colorEntity, rank );
      }
    }
  }
}


static void
EG_getBodyUnits(const egObject *model, int nBody, const egObject **bodies,
                const char **units)
{
  int          i, stat, aType, aLen;
  const int    *ints;
  const double *reals;
  const char   *str;

  *units = NULL;
  if (model != NULL) {
    stat = EG_attributeRet(model, ".lengthUnits", &aType, &aLen, &ints,
                           &reals, &str);
    if (stat == EGADS_SUCCESS)
      if (aType == ATTRSTRING) {
        *units = str;
        return;
      }
  }

  for (i = 0; i < nBody; i++) {
    stat = EG_attributeRet(bodies[i], ".lengthUnits", &aType, &aLen, &ints,
                           &reals, &str);
    if (stat == EGADS_SUCCESS)
      if (aType == ATTRSTRING)
        if (*units == NULL) {
          *units = str;
        } else {
          if (strcmp(*units, str) != 0) {
            printf(" EGADS Warning: Inconsistent Units (EG_saveModel)!\n");
            *units = NULL;
            return;
          }
        }
  }

}


int
EG_saveModel(const egObject *model, const char *name)
{
  int            i, j, n, len, outLevel, ibody, nbody, stat;
  double         scale = 1.0;
  TopoDS_Shape   wshape;
  const egObject *mdl, **objs;
  const char     *units, *wunits;
  FILE           *fp;

  if  (model == NULL)               return EGADS_NULLOBJ;
  if  (model->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if ((model->oclass != MODEL) &&
      (model->oclass != BODY))      return EGADS_NOTMODEL;
  outLevel = EG_outLevel(model);

  if (name == NULL) {
    if (outLevel > 0)
      printf(" EGADS Warning: NULL Filename (EG_saveModel)!\n");
    return EGADS_NONAME;
  }

  /* does file exist? */

  fp = fopen(name, "r");
  if (fp != NULL) {
    if (outLevel > 0)
      printf(" EGADS Warning: File %s Exists (EG_saveModel)!\n", name);
    fclose(fp);
    return EGADS_EXISTS;
  }

  /* find extension */

  len = strlen(name);
  for (i = len-1; i > 0; i--)
    if (name[i] == '.') break;
  if (i == 0) {
    if (outLevel > 0)
      printf(" EGADS Warning: No Extension in %s (EG_saveModel)!\n", name);
    return EGADS_NODATA;
  }

  if (model->oclass == MODEL) {
    egadsModel *mshape = (egadsModel *) model->blind;
    wshape             = mshape->shape;
    nbody              = mshape->nbody;
    objs               = (const egObject **) mshape->bodies;
    mdl                = model;
  } else {
    egadsBody *mshape  = (egadsBody *) model->blind;
    wshape             = mshape->shape;
    nbody              = 1;
    objs               = &model;
    mdl                = NULL;
  }

  if ((strcasecmp(&name[i],".step") == 0) ||
      (strcasecmp(&name[i],".stp") == 0)) {

    /* STEP files */

    EG_getBodyUnits(mdl, nbody, objs, &units);
    if (units == NULL) {
      wunits = "MM";
    } else {
      EG_importScale("STEP Writer", units, &scale, &wunits);
    }

    STEPControl_Writer aWriter;
    Interface_Static::SetCVal("write.step.unit", wunits);
    const Handle(XSControl_WorkSession)& theSession = aWriter.WS();
    const Handle(XSControl_TransferWriter)& aTransferWriter =
                                                  theSession->TransferWriter();
    const Handle(Transfer_FinderProcess) FP = aTransferWriter->FinderProcess();
    APIHeaderSection_MakeHeader makeHeader(aWriter.Model());
    Handle(TCollection_HAsciiString) headerAuthor = new
                                          TCollection_HAsciiString("ESP EGADS");
    makeHeader.SetAuthorValue(1, headerAuthor);

    TopExp_Explorer Exp;
    const STEPControl_StepModelType aVal = STEPControl_AsIs;
    for (Exp.Init(wshape, TopAbs_WIRE,  TopAbs_FACE);
         Exp.More(); Exp.Next()) aWriter.Transfer(Exp.Current(), aVal);
    for (Exp.Init(wshape, TopAbs_FACE,  TopAbs_SHELL);
         Exp.More(); Exp.Next()) aWriter.Transfer(Exp.Current(), aVal);
    for (Exp.Init(wshape, TopAbs_SHELL, TopAbs_SOLID);
         Exp.More(); Exp.Next()) aWriter.Transfer(Exp.Current(), aVal);
    for (Exp.Init(wshape, TopAbs_SOLID);
         Exp.More(); Exp.Next()) aWriter.Transfer(Exp.Current(), aVal);

    // transfer Name attributes to StepRepr_RepresentationItem
    EG_setSTEPprops(theSession, nbody, objs, FP);

    if (!aWriter.Write(name)) {
      printf(" EGADS Warning: STEP Write Error (EG_saveModel)!\n");
      return EGADS_WRITERR;
    }

  } else if ((strcasecmp(&name[i],".iges") == 0) ||
             (strcasecmp(&name[i],".igs") == 0)) {

    /* IGES files */

    egadsShapeData shapeData;
    EG_getShapeData(nbody, objs, shapeData);

    EG_getBodyUnits(mdl, nbody, objs, &units);
    if (units == NULL) {
      wunits = "MM";
    } else {
      EG_importScale("IGES Writer", units, &scale, &wunits);
      if (scale != 1.0) {
        scale = 1.0/scale;
        gp_Trsf form = gp_Trsf();
        form.SetValues(scale, 0.0,   0.0,   0.0,
                       0.0,   scale, 0.0,   0.0,
                       0.0,   0.0,   scale, 0.0);
        BRepBuilderAPI_Transform xForm(wshape, form, Standard_True);
        if (!xForm.IsDone()) {
          printf(" EGADS Warning: Can't scale Object (EG_saveModel)!\n");
        } else {
          wshape = shapeData.Update(wshape, xForm);
        }
      }
    }

    try {
      IGESControl_Controller::Init();
      Standard_Integer modecr = 1;    // BRep export
      IGESControl_Writer iWrite(wunits, modecr);
      TopExp_Explorer Exp;
      for (Exp.Init(wshape, TopAbs_WIRE,  TopAbs_FACE);
           Exp.More(); Exp.Next()) iWrite.AddShape(Exp.Current());
      for (Exp.Init(wshape, TopAbs_FACE,  TopAbs_SHELL);
           Exp.More(); Exp.Next()) iWrite.AddShape(Exp.Current());
      for (Exp.Init(wshape, TopAbs_SHELL, TopAbs_SOLID);
           Exp.More(); Exp.Next()) iWrite.AddShape(Exp.Current());
      for (Exp.Init(wshape, TopAbs_SOLID);
           Exp.More(); Exp.Next()) iWrite.AddShape(Exp.Current());

      // transfer Name/Color attributes to IGES
      EG_setIGESprop(shapeData, iWrite);

      iWrite.ComputeModel();
      if (!iWrite.Write(name)) {
        printf(" EGADS Error: IGES Write Error (EG_saveModel)!\n");
        return EGADS_WRITERR;
      }
    }
    catch (...)
    {
      printf(" EGADS Error: Internal IGES Write Error (EG_saveModel)!\n");
      return EGADS_WRITERR;
    }

  } else if ((strcasecmp(&name[i],".brep") == 0) ||
             (strcasecmp(&name[i],".egads") == 0)) {

    /* Native OCC file or our filetype */

#if CASVER < 760
    if (!BRepTools::Write(wshape, name)) {
#else
    if (!BRepTools::Write(wshape, name, Standard_False, Standard_False,
                          TopTools_FormatVersion_VERSION_2)) {
#endif
      printf(" EGADS Warning: OCC Write Error (EG_saveModel)!\n");
      return EGADS_WRITERR;
    }
    if (strcasecmp(&name[i],".brep") == 0) return EGADS_SUCCESS;

    /* append the attributes -- output in the read order */

    FILE *fp = fopen(name, "a");
    if (fp == NULL) {
      printf(" EGADS Warning: EGADS Open Error (EG_saveModel)!\n");
      return EGADS_WRITERR;
    }
    fprintf(fp, "\n##EGADS HEADER FILE-REV 1 ##\n");
    /* write model attributes */
    if (model->oclass == MODEL) {
      EG_writeAttrs(model, fp);
    } else {
      fprintf(fp, "0\n");
    }
    TopExp_Explorer Exp;
    for (Exp.Init(wshape, TopAbs_WIRE,  TopAbs_FACE); Exp.More(); Exp.Next()) {
      TopoDS_Shape shape = Exp.Current();
      for (i = 0; i < nbody; i++) {
        const egObject *obj   = objs[i];
        egadsBody      *pbody = (egadsBody *) obj->blind;
        if (shape.IsSame(pbody->shape)) {
          EG_writeAttrs(obj, fp);
          break;
        }
      }
    }
    for (Exp.Init(wshape, TopAbs_FACE,  TopAbs_SHELL); Exp.More(); Exp.Next()) {
      TopoDS_Shape shape = Exp.Current();
      for (i = 0; i < nbody; i++) {
        const egObject *obj   = objs[i];
        egadsBody      *pbody = (egadsBody *) obj->blind;
        if (shape.IsSame(pbody->shape)) {
          EG_writeAttrs(obj, fp);
          break;
        }
      }
    }
    for (Exp.Init(wshape, TopAbs_SHELL, TopAbs_SOLID); Exp.More(); Exp.Next()) {
      TopoDS_Shape shape = Exp.Current();
      for (i = 0; i < nbody; i++) {
        const egObject *obj   = objs[i];
        egadsBody      *pbody = (egadsBody *) obj->blind;
        if (shape.IsSame(pbody->shape)) {
          EG_writeAttrs(obj, fp);
          break;
        }
      }
    }
    for (Exp.Init(wshape, TopAbs_SOLID); Exp.More(); Exp.Next()) {
      TopoDS_Shape shape = Exp.Current();
      for (i = 0; i < nbody; i++) {
        const egObject *obj   = objs[i];
        egadsBody      *pbody = (egadsBody *) obj->blind;
        if (shape.IsSame(pbody->shape)) {
          EG_writeAttrs(obj, fp);
          break;
        }
      }
    }
    /* add non-Body types */
    if ((model->oclass == MODEL) && (model->mtype > nbody)) {
      fprintf(fp, "##EGADS HEADER FILE-EXT 1 ##\n");
      fprintf(fp, "%hd\n", model->mtype);
      for (j = nbody; j < model->mtype; j++) {
        const egObject *obj = objs[j];
        TopoDS_Shape bshape;
        if (obj->oclass == TESSELLATION) {
          egTessel  *btess = (egTessel *) obj->blind;
          egObject  *src   = btess->src;
          if (src->oclass == EBODY) {
            for (i = nbody; i < model->mtype; i++) {
              const egObject *bod = objs[i];
              if (bod == src) break;
            }
            if (i == model->mtype) {
              printf(" EGADS Internal: Ancillary tess %d -- cannot find EBody!\n",
                     j+1);
              break;
            }
            fprintf(fp, "%d %d\n", obj->oclass, i+1);
            stat = EG_writeTess(obj, fp);
            if (stat != EGADS_SUCCESS) {
              printf(" EGADS Error: Ancillary tessellation %d -- status = %d\n",
                     j+1, stat);
              break;
            }
            continue;
          }
          egadsBody *pbody = (egadsBody *) src->blind;
          bshape           = pbody->shape;
        } else {
          egEBody   *ebody = (egEBody *) obj->blind;
          egObject  *ref   = ebody->ref;
          egadsBody *pbody = (egadsBody *) ref->blind;
          bshape           = pbody->shape;
        }
        n = ibody = 0;
        for (Exp.Init(wshape, TopAbs_WIRE,  TopAbs_FACE); Exp.More(); Exp.Next()) {
          TopoDS_Shape shape = Exp.Current();
          for (i = 0; i < nbody; i++) {
            const egObject *obj   = objs[i];
            egadsBody      *pbody = (egadsBody *) obj->blind;
            if (shape.IsSame(pbody->shape)) {
              n++;
              if (shape.IsSame(bshape)) ibody = n;
              break;
            }
          }
        }
        for (Exp.Init(wshape, TopAbs_FACE,  TopAbs_SHELL); Exp.More(); Exp.Next()) {
          TopoDS_Shape shape = Exp.Current();
          for (i = 0; i < nbody; i++) {
            const egObject *obj   = objs[i];
            egadsBody      *pbody = (egadsBody *) obj->blind;
            if (shape.IsSame(pbody->shape)) {
              n++;
              if (shape.IsSame(bshape)) ibody = n;
              break;
            }
          }
        }
        for (Exp.Init(wshape, TopAbs_SHELL, TopAbs_SOLID); Exp.More(); Exp.Next()) {
          TopoDS_Shape shape = Exp.Current();
          for (i = 0; i < nbody; i++) {
            const egObject *obj   = objs[i];
            egadsBody      *pbody = (egadsBody *) obj->blind;
            if (shape.IsSame(pbody->shape)) {
              n++;
              if (shape.IsSame(bshape)) ibody = n;
              break;
            }
          }
        }
        for (Exp.Init(wshape, TopAbs_SOLID); Exp.More(); Exp.Next()) {
          TopoDS_Shape shape = Exp.Current();
          for (i = 0; i < nbody; i++) {
            const egObject *obj   = objs[i];
            egadsBody      *pbody = (egadsBody *) obj->blind;
            if (shape.IsSame(pbody->shape)) {
              n++;
              if (shape.IsSame(bshape)) ibody = n;
              break;
            }
          }
        }
        if ((n != nbody) || (ibody == 0)) {
          printf(" EGADS Internal: Ancillary objects -- n = %d [%d], ib = %d\n",
                 n, nbody, ibody);
          break;
        }
        fprintf(fp, "%d %d\n", obj->oclass, ibody);
        if (obj->oclass == TESSELLATION) {
          stat = EG_writeTess(obj, fp);
        } else {
          stat = EG_writeEBody(obj, fp);
        }
        if (stat != EGADS_SUCCESS) {
          printf(" EGADS Error: Ancillary objects %d -- status = %d\n",
                 j+1, stat);
          break;
        }
      }
    }
    fclose(fp);

  } else {
    if (outLevel > 0)
      printf(" EGADS Warning: Extension in %s Not Supported (EG_saveModel)!\n",
             name);
    return EGADS_NODATA;
  }

  return EGADS_SUCCESS;
}


int
EG_saveTess(egObject *tess, const char *name)
{
  int          status, len,   outLevel, stat;
  int          npts,   ntri,  iedge,    iface;
  int          nnode,  nedge, nface,    nattr = 0;
  const double *pxyz  = NULL, *puv    = NULL, *pt    = NULL;
  const int    *ptype = NULL, *pindex = NULL, *ptris = NULL, *ptric = NULL;
  egObject     *body;
  egAttrs      *attrs;

  if (tess == NULL)                 return EGADS_NULLOBJ;
  if (tess->magicnumber != MAGIC)   return EGADS_NOTOBJ;
  if (tess->oclass != TESSELLATION) return EGADS_NOTTESS;
  outLevel = EG_outLevel(tess);

  /* does file exist? */

  FILE *fp = fopen(name, "rb");
  if (fp != NULL) {
    if (outLevel > 0)
      printf(" EGADS Warning: File %s Exists (EG_saveTess)!\n", name);
    status = EGADS_EXISTS;
    goto cleanup;
  }

  // get the body from tessellation
  status = EG_statusTessBody(tess, &body, &stat, &npts);
  if (status != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Warning: Tessellation is Open (EG_saveTess)!\n");
    status = EGADS_TESSTATE;
    goto cleanup;
  }

  status = EG_getBodyTopos(body, NULL, NODE, &nnode, NULL);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_getBodyTopos(body, NULL, EDGE, &nedge, NULL);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_getBodyTopos(body, NULL, FACE, &nface, NULL);
  if (status != EGADS_SUCCESS) goto cleanup;

  fp = fopen(name, "wb");
  if (fp == NULL) {
    printf(" EGADS Warning: File %s Open Error (EG_saveTess)!\n", name);
    status = EGADS_WRITERR;
    goto cleanup;
  }

  // write the number of nodes, edges, and faces
  fwrite(&nnode, sizeof(int), 1, fp);
  fwrite(&nedge, sizeof(int), 1, fp);
  fwrite(&nface, sizeof(int), 1, fp);

  // write out the edge tessellation
  for (iedge = 0; iedge < nedge; iedge++) {

    status = EG_getTessEdge(tess, iedge+1, &len, &pxyz, &pt);
    if (status != EGADS_SUCCESS) goto cleanup;

    fwrite(&len, sizeof(int),        1, fp);
    if (len == 0) continue;
    fwrite(pxyz, sizeof(double), 3*len, fp);
    fwrite(pt,   sizeof(double),   len, fp);
    if (ferror(fp)) {
      printf ("EGADS Warning: File %s Write Error (EG_saveTess)!\n", name);
      status = EGADS_WRITERR;
      goto cleanup;
    }
  }

  // write out face tessellations
  for (iface = 0; iface < nface; iface++) {
    EG_getTessFace(tess, iface+1, &len, &pxyz, &puv, &ptype, &pindex,
                   &ntri, &ptris, &ptric);

    fwrite(&len,   sizeof(int),         1, fp);
    fwrite(&ntri,  sizeof(int),         1, fp);
    if ((len == 0) || (ntri == 0)) continue;
    fwrite(pxyz,   sizeof(double),  3*len, fp);
    fwrite(puv,    sizeof(double),  2*len, fp);
    fwrite(ptype,  sizeof(int),       len, fp);
    fwrite(pindex, sizeof(int),       len, fp);
    fwrite(ptris,  sizeof(int),    3*ntri, fp);
    fwrite(ptric,  sizeof(int),    3*ntri, fp);
    if (ferror(fp)) {
      printf ("EGADS Warning: File %s Write Error (EG_saveTess)!\n", name);
      status = EGADS_WRITERR;
      goto cleanup;
    }
  }

  // write out the tessellation attributes
  attrs = (egAttrs *) tess->attrs;
  if (attrs != NULL) nattr = attrs->nattrs;
  fprintf(fp, "%d\n", nattr);
  if (nattr != 0) EG_writeAttr(attrs, fp);

  status = EGADS_SUCCESS;

cleanup:
  fclose(fp);

  return status;
}


int
EG_loadTess(egObject *body, const char *name, egObject **tess)
{
  int      i, ir, status, nnode, nedge, nface, n[3], len, ntri, nattr;
  int      *ptype, *pindex, *tris, *tric;
  double   *xyz, *param;
  egTessel *btess;
  egObject *obj;
  FILE     *fp;

  *tess  = NULL;
  status = EG_getBodyTopos(body, NULL, NODE, &nnode, NULL);
  if (status != EGADS_SUCCESS) return status;
  if (body->oclass == EBODY) {
    status = EG_getBodyTopos(body, NULL, EEDGE, &nedge, NULL);
    if (status != EGADS_SUCCESS) return status;
    status = EG_getBodyTopos(body, NULL, EFACE, &nface, NULL);
    if (status != EGADS_SUCCESS) return status;
  } else {
    status = EG_getBodyTopos(body, NULL, EDGE, &nedge, NULL);
    if (status != EGADS_SUCCESS) return status;
    status = EG_getBodyTopos(body, NULL, FACE, &nface, NULL);
    if (status != EGADS_SUCCESS) return status;
  }

  fp = fopen(name, "rb");
  if (fp == NULL) {
    printf(" EGADS Error: File %s Does not Exist (EG_loadTess)!\n", name);
    return EGADS_EXISTS;
  }

  ir = fread(n, sizeof(int), 3, fp);
  if (ir != 3) {
    printf(" EGADS Error: Header with only %d words (EG_loadTess)!\n", ir);
    fclose(fp);
    return EGADS_INDEXERR;
  }
  if ((nnode != n[0]) || (nedge != n[1]) || (nface != n[2])) {
    printf(" EGADS Error: Count mismatch %d %d  %d %d  %d %d (EG_loadTess)!\n",
           nnode, n[0], nedge, n[1], nface, n[2]);
    fclose(fp);
    return EGADS_INDEXERR;
  }

  /* initialize the Tessellation Object */
  status = EG_initTessBody(body, tess);
  if (status != EGADS_SUCCESS) {
    fclose(fp);
    return status;
  }

  /* do the Edges */
  for (i = 0; i < nedge; i++) {
    len = 0;
    fread(&len, sizeof(int), 1, fp);
    if (len == 0) continue;
    xyz   = (double *) malloc(3*len*sizeof(double));
    param = (double *) malloc(  len*sizeof(double));
    if ((xyz == NULL) || (param == NULL)) {
      printf(" EGADS Error: malloc on Edge %d -- len = %d (EG_loadTess)!\n",
               i+1, len);
      if (xyz   != NULL) free(xyz);
      if (param != NULL) free(param);
      EG_deleteObject(*tess);
      *tess = NULL;
      fclose(fp);
      return EGADS_MALLOC;
    }
    fread(xyz,   sizeof(double), 3*len, fp);
    fread(param, sizeof(double),   len, fp);
    status = EG_setTessEdge(*tess, i+1, len, xyz, param);
    free(xyz);
    free(param);
    if (status != EGADS_SUCCESS) {
      printf(" EGADS Error: EG_setTessEdge %d = %d (EG_loadTess)!\n",
             i+1, status);
      EG_deleteObject(*tess);
      *tess = NULL;
      fclose(fp);
      return status;
    }
  }

  /* do the Faces */
  for (i = 0; i < nface; i++) {
    len = ntri = 0;
    fread(&len,  sizeof(int), 1, fp);
    fread(&ntri, sizeof(int), 1, fp);
    if ((len == 0) || (ntri == 0)) {
      btess = (egTessel *) (*tess)->blind;
      btess->nFace = 0;
      EG_free(btess->tess2d);
      btess->tess2d = NULL;
      continue;
    }
    xyz    = (double *) malloc(3*len*sizeof(double));
    param  = (double *) malloc(2*len*sizeof(double));
    ptype  = (int *)    malloc(  len* sizeof(int));
    pindex = (int *)    malloc(  len* sizeof(int));
    tris   = (int *)    malloc(3*ntri*sizeof(int));
    tric   = (int *)    malloc(3*ntri*sizeof(int));
    if ((xyz    == NULL) || (param == NULL) || (ptype == NULL) ||
        (pindex == NULL) || (tris  == NULL) || (tric  == NULL)) {
      printf(" EGADS Error: malloc on Face %d -- lens = %d %d (EG_loadTess)!\n",
             i+1, len, ntri);
      if (xyz    != NULL) free(xyz);
      if (param  != NULL) free(param);
      if (ptype  != NULL) free(ptype);
      if (pindex != NULL) free(pindex);
      if (tris   != NULL) free(tris);
      if (tric   != NULL) free(tric);
      EG_deleteObject(*tess);
      *tess = NULL;
      fclose(fp);
      return EGADS_MALLOC;
    }
    fread(xyz,    sizeof(double), 3*len,  fp);
    fread(param,  sizeof(double), 2*len,  fp);
    fread(ptype,  sizeof(int),      len,  fp);
    fread(pindex, sizeof(int),      len,  fp);
    fread(tris,   sizeof(int),    3*ntri, fp);
    fread(tric,   sizeof(int),    3*ntri, fp);
    status = EG_setTessFace(*tess, i+1, len, xyz, param, ntri, tris);
    free(xyz);
    free(param);
    free(ptype);
    free(pindex);
    free(tris);
    free(tric);
    if (status != EGADS_SUCCESS) {
      printf(" EGADS Error: EG_setTessFace %d = %d (EG_loadTess)!\n",
             i+1, status);
      EG_deleteObject(*tess);
      *tess = NULL;
      fclose(fp);
      return status;
    }
  }

  /* close up the open tessellation */
  status = EG_statusTessBody(*tess, &obj, &i, &len);
  if (status != EGADS_SUCCESS) {
    printf(" EGADS Error: EG_statusTessBody = %d (EG_loadTess)!\n", status);
    EG_deleteObject(*tess);
    *tess = NULL;
    fclose(fp);
    return status;
  }
  if (i != 1) {
    printf(" EGADS Error: New Tessellation Object is Open (EG_loadTess)!\n");
    EG_deleteObject(*tess);
    *tess = NULL;
    fclose(fp);
    return EGADS_TESSTATE;
  }

  /* attach the attributes */
  fscanf(fp, "%d\n", &nattr);
  if (nattr != 0) EG_readAttrs(*tess, nattr, fp);
  fclose(fp);

  return EGADS_SUCCESS;
}

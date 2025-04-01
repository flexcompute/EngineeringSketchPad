#ifndef EGADSCLASSES_H
#define EGADSCLASSES_H
/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             C++/OpenCASCADE Object Header
 *
 *      Copyright 2011-2025, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include "egadsOCC.h"
#include "Surreal/SurrealS.h"

class egadsPCurve
{
public:
  Handle(Geom2d_Curve) handle;
  egObject             *ref;
  int                  topFlg;
  int                  dataLen;
  int                  *header;
  double               *data;
  SurrealS<1>          *data_dot;
  double               trange[2];
};


class egadsCurve
{
public:
  Handle(Geom_Curve) handle;
  egObject           *ref;
  int                topFlg;
  int                dataLen;
  int                *header;
  double             *data;
  SurrealS<1>        *data_dot;
  double             trange[2];
};


class egadsSurface
{
public:
  Handle(Geom_Surface) handle;
  egObject             *ref;
  int                  topFlg;
  int                  dataLen;
  int                  *header;
  double               *data;
  SurrealS<1>          *data_dot;
  double               urange[2];
  double               vrange[2];
};


class egadsBox
{
public:
  int    filled;
  double box[6];
};


class egadsNode
{
public:
  TopoDS_Vertex node;
  double        xyz[3];
  egadsBox      bbox;
  int           filled;
  SurrealS<1>   xyz_dot[3];
};


class egadsEdge
{
public:
  TopoDS_Edge edge;
  egObject    *curve;                   // curve object
  egObject    *nodes[2];                // pointer to ego nodes
  int         topFlg;
  double      trange[2];
  egadsBox    bbox;
  int         filled;
  SurrealS<1> trange_dot[2];
};


class egadsLoop
{
public:
  TopoDS_Wire loop;
  egObject    *surface;                 // associated non-planar surface
                                        // will have pcurves after edges (nonNULL)
  int         nedges;                   // number of edges
  egObject    **edges;                  // edge objects (*2 if surface is nonNULL)
  int         *senses;                  // sense for each edge
  int         topFlg;
  egadsBox    bbox;
};


class egadsFace
{
public:
  TopoDS_Face face;
  egObject    *surface;                 // surface object
  int         nloops;                   // number of loops
  egObject    **loops;                  // loop objects
  int         *senses;                  // outer/inner for each loop
  int         topFlg;
  double      urange[2];
  double      vrange[2];
  egadsBox    bbox;
};


class egadsShell
{
public:
  TopoDS_Shell shell;
  int          nfaces;                  // number of faces
  egObject    **faces;                  // face objects
  int         topFlg;
  egadsBox    bbox;
};


class egadsMap
{
public:
  TopTools_IndexedMapOfShape map;
  egObject                   **objs;    // vector of egos w/ map
};


class egadsBody
{
public:
  TopoDS_Shape shape;                   // OCC topology
  egadsMap     nodes;
  egadsMap     edges;
  egadsMap     loops;
  egadsMap     faces;
  egadsMap     shells;
  int          *senses;                 // shell outer/inner (solids)
  egadsBox     bbox;
  int          massFill;
  double       massProp[14];
};


class egadsModel
{
public:
  TopoDS_Shape shape;                   // OCC topology
  int          nobjs;                   // number of total egObjects
  int          nbody;                   // number of bodies
  egObject     **bodies;                // vector of pointers to egObjects
  egadsBox     bbox;
};


// Used to track labels from STEP/IGES readers
class egadsShapeData
{
  class Label
  {
  public:
    Label(const char* shapeName) : shapeName(EG_strdup(shapeName)) {}
    Label(const Label& label) : shapeName(EG_strdup(label.shapeName)) {}
    ~Label() { EG_free(shapeName); }

    char  *shapeName;
  };
  class BodyColor
  {
  public:
    BodyColor(const Quantity_Color* Edge,
              const Quantity_Color* Face) : Edge(Edge), Face(Face) {}
    BodyColor(const BodyColor& bodycolor) :
      Edge( bodycolor.Edge != NULL ? new Quantity_Color(*bodycolor.Edge) : NULL ),
      Face( bodycolor.Face != NULL ? new Quantity_Color(*bodycolor.Face) : NULL ) {}
    BodyColor& operator=(const BodyColor& bodycolor) {
      Edge = bodycolor.Edge != NULL ? new Quantity_Color(*bodycolor.Edge) : NULL;
      Face = bodycolor.Face != NULL ? new Quantity_Color(*bodycolor.Face) : NULL;
      return *this;
    }
    ~BodyColor() { delete Edge; delete Face; }

    const Quantity_Color *Edge;
    const Quantity_Color *Face;
  };
  typedef NCollection_IndexedDataMap<TopoDS_Shape, Label, TopTools_ShapeMapHasher> Label_IndexedDataMap;
  typedef NCollection_IndexedDataMap<TopoDS_Shape, Quantity_Color, TopTools_ShapeMapHasher> Color_IndexedDataMap;
  typedef NCollection_IndexedDataMap<TopoDS_Shape, BodyColor, TopTools_ShapeMapHasher> BodyColor_IndexedDataMap;
public:

  egadsShapeData() {}

  void Add(const TopoDS_Shape& shape, const char* shapeName) { labels.Add(shape, Label(shapeName)); }
  void Add(const TopoDS_Shape& shape, const Quantity_Color& color) { colors.Add(shape, color); }
  void Add(const TopoDS_Shape& shape, const Quantity_Color* Edge,
                                      const Quantity_Color* Face) { bodycolors.Add(shape, BodyColor(Edge,Face)); }

  Standard_Integer LabelExtent() const { return labels.Extent(); }
  Standard_Integer LabelFindIndex(const TopoDS_Shape& shape) const { return labels.FindIndex(shape); }
  const TopoDS_Shape& LabelFindKey(Standard_Integer i) const { return labels.FindKey(i); }
  const Label& label(Standard_Integer i) const { return labels(i); }

  Standard_Integer ColorExtent() const { return colors.Extent(); }
  Standard_Integer ColorFindIndex(const TopoDS_Shape& shape) const { return colors.FindIndex(shape); }
  const TopoDS_Shape& ColorFindKey(Standard_Integer i) const { return colors.FindKey(i); }
  const Quantity_Color& color(Standard_Integer i) const { return colors(i); }

  Standard_Integer BodyColorExtent() const { return bodycolors.Extent(); }
  Standard_Integer BodyColorFindIndex(const TopoDS_Shape& shape) const { return bodycolors.FindIndex(shape); }
  const TopoDS_Shape& BodyColorFindKey(Standard_Integer i) const { return bodycolors.FindKey(i); }
  const Quantity_Color* EdgeColor(Standard_Integer i) const { return bodycolors(i).Edge; }
  const Quantity_Color* FaceColor(Standard_Integer i) const { return bodycolors(i).Face; }

  TopoDS_Shape Update(const TopoDS_Shape& oldShape, const TopoDS_Shape& newShape);
  TopoDS_Shape Update(const TopoDS_Shape& oldShape, BRepBuilderAPI_ModifyShape& xForm);
  TopoDS_Shape Update(const TopoDS_Shape& oldShape, const Handle(BRepTools_ReShape)& reShape);
  TopoDS_Shape Update(const TopoDS_Shape& oldShape, const TopoDS_Shape& newShape, const Handle(BRepTools_History)& history);

  void AddShapeName(egadsMap& shapes);
  void AddShapeColor(egadsMap& shapes);
protected:
  Label_IndexedDataMap labels;
  Color_IndexedDataMap colors;
  BodyColor_IndexedDataMap bodycolors;
};




#endif

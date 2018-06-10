/*
CREATED : Jan 31, 2013
MODIFIED:
AUTHOR  : Yu-Chung Hsiao
EMAIL   : project.caplet@gmail.com

This file is part of CAPLET.

CAPLET is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

CAPLET is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the Lesser GNU General Public License
along with CAPLET.  If not, see <http://www.gnu.org/licenses/>.
*/


#ifndef GDSGEOMETRY_H
#define GDSGEOMETRY_H

#include <list>
#include <cmath>
#include <map>
#include <vector>
#include <stdexcept>
#include <functional>
#include <sstream>


//****
//*
//* This file defines basic geometry data structures used by
//* GeoLoader and PanelRenderer
//*
//****

//****
//*
//* Exceptions
//*
//*

class GeometryNotManhattanError : public std::logic_error
{
public:
    explicit GeometryNotManhattanError()
        : logic_error(""){}
    virtual ~GeometryNotManhattanError() throw() {}
    virtual const char* what() const throw(){
        return "Geometries are not Manhattan.";
    }
};

class ShapeTransformationError : public std::logic_error
{
public:
    explicit ShapeTransformationError (const std::string msg)
        : logic_error(msg), m_msg(msg) {}
    virtual ~ShapeTransformationError() throw() {}
    virtual const char* what() const throw(){
        return m_msg.c_str();
    }
private:
    std::string m_msg;
};


class ConductorLayerNotCompatibleError : public std::logic_error{
public:
    explicit ConductorLayerNotCompatibleError (
            const unsigned nMetal1,
            const unsigned nVia1,
            const unsigned nMetal2,
            const unsigned nVia2) : logic_error(""){

        std::stringstream ss1;
        ss1 << "(" << nMetal1 << ", " << nVia1 << ")";
        m_layerInfo1 = ss1.str();
        std::stringstream ss2;
        ss2 << "(" << nMetal2 << ", " << nVia2 << ")";
        m_layerInfo2 = ss2.str();
    }
    virtual ~ConductorLayerNotCompatibleError() throw() {}
    virtual const char* what() const throw(){
        return ("Conductor layer not compatible: " + m_layerInfo1 + " != " + m_layerInfo2).c_str();
    }
private:
    std::string m_layerInfo1;
    std::string m_layerInfo2;
};


//****
//*
//* Geometries
//*
//*

//**
//* Dir
enum Dir {  X=1, Y=2, Z=3, XP=1, XM=-1, YP=2, YM=-2, ZP=3, ZM=-3, FLAT=0 };
const int nDir = 6;
//**
//* Point
class Point{
public:
    int x;
    int y;
    int z;
    Dir dir;
    int len; //* only used in a polygon


    Point();
    Point(int xx, int yy);
    Point(int xx, int yy, int zz);

    //@author:Dragon
    //add a method to Point class: setCoordinate();
    void setCoordinate(int xx, int yy, int zz );

    Point operator-(const Point& other);
    bool  operator==(const Point& other);
    inline int vecLen2(){ return (x!=0)? int(std::abs(float(x))) : int(std::abs(float(y))); } //* 2D Manhattan geometry only
};

//****
//*
//* Polygon and its collection
//*
//*

//**
//* Polygon
//typedef std::list<Point> Polygon;

typedef std::list<Point> PointList;

class Polygon : public PointList{
public:
    explicit Polygon( const allocator_type &allo=allocator_type());
    explicit Polygon(
             size_type              n,
             const Point            &value = Point(),
             const allocator_type   &allo=allocator_type());
    Polygon( iterator               first,
             iterator               last,
             const allocator_type   &allo=allocator_type());
    Polygon( const Polygon& poly);

    bool isManhattan() const;
};

//**
//* PolygonList
typedef std::list<Polygon> PolygonList;

//**
//* LayeredPolygonList
typedef std::vector<PolygonList> LayeredPolygonList;


//**
//* Rectangle
class Rectangle{
public:
    Dir     normal;
    int     x1;
    int     x2;
    int     y1;
    int     y2;
    int     z1;
    int     z2;

    //**
    //* Default constructor
    Rectangle();
    Rectangle(Dir normal, int x1, int x2, int y1, int y2, int z1, int z2);

    //@author:Dragon
    Rectangle(  int x1, int x2, int y1, int y2, int z1, int z2 );

    //**
    //* Constructor
    //* - Transform a 2D (+z-dir) 4-point polygon to a rectangle
    Rectangle( const Polygon &poly) throw (ShapeTransformationError);

    //**
    //* Copy constructor
    Rectangle( const Rectangle &rect );

    //**
    //* area
    //* - return area of the rectangle
    double  area() const;

    //**
    //* geometry difference
    Rectangle operator- (const Rectangle &rect ) const;

    //**
    //* isOverlapping
    //* - return true if this overlaps rect considering only x and y
    bool    isOverlapping  ( const Rectangle &rect ) const;
    bool    isOverlapping3d( const Rectangle &rect ) const;

    //@author:Dragon
    //to check if two cuboid is insecting
    bool isIntersecting3d ( const Rectangle &rect ) const;

    //@author:Dragon
    //to check if a rect ( 2d or 3d) is within a rect (2d or 3d)
    bool isIncludedIn( const Rectangle &rect ) const;

    //**
    //* hasCornerInside
    //* - return true if this has a corner in rect
    bool    hasCornerInside ( const Rectangle &rect ) const;

    //**
    //* print
    //* - print direction, area, and coordinates
    void print() const;

//@author:Dragon
//This fuction is to return the max length of xLen,yLen,and zLen
    int size() const;

};


//****
//*
//* Rectangle Collection
//*
//*

//**
//* RectangleList
class RectangleList : public std::list<Rectangle>
{
public:
    explicit RectangleList(
            const allocator_type    &allo = allocator_type());
    explicit RectangleList(
            size_type               n,
            const Rectangle         &value = Rectangle(),
            const allocator_type    &allo = allocator_type());
    RectangleList(
            iterator                first,
            iterator                last,
            const allocator_type    &allo = allocator_type());
    RectangleList( const RectangleList& rectList);

    //**
    //* merge
    //* - self merge
    void merge();

    //**
    //* decompose
    //* - make this disjoint rect set
    void decompose();
};

//**
//* LayeredRectangleList
typedef std::vector<RectangleList> LayeredRectangleList;

//**
//* ConnectedRectangleList
typedef std::list<RectangleList> ConnectedRectangleList;

//**
//* LayeredConnectedRectangleList
typedef std::vector<ConnectedRectangleList> LayeredConnectedRectangleList;


//**
//* RectangleMap
//* - for incremental sorting
typedef std::multimap<double, Rectangle, std::greater<double> > RectangleMap;



//@author:Dragon
//here we place the notification of the GaussianSurface
class GaussianSurface;
class GaussianSurfaceList:public std::list<GaussianSurface>{
public:
    //@author:Dragon
    //the effective area of the virtual Gaussian Surface of the conductor
    double VGSArea;
    //@author:Dragon
    //the surface integral of the PDF on the gaussian surface which can be approximated by sampling a point on gaussian surface and compute its PDF value
    //then multiplied by VGSArea
    double PDFIntegralOnVGS;

    //@author:Dragon
    //the sum of block Gaussian Surface Area
    double SumOfBGSArea;
    //Nt: total number of the points generated
    //Ng: the number points which locate exactly on VGS uniformly
    long Nt,Ng;

    //the minimum extension Distance of the gaussian surface, which is used to generate the up limit of the PDF of gaussian surface
    double minExtensionDis;

    GaussianSurfaceList(){
        VGSArea=0;
        //This variable is H in the formula
        PDFIntegralOnVGS=0;
        SumOfBGSArea=0;
        Nt=0;
        Ng=0;
        minExtensionDis=0;
    }

};



//@author:Dragon
/*
SubConductor is the kid Class of Rectangule. It is used to express the rectangulars which are included in a conductor, and it also contains the pointer of 
that conductor (fatherConductor), this is an important property----------This is wrong!

*/
class SubConductor: public Rectangle{
public:
    int fatherConductorID;

    //@author:Dragon
    //用该变量判断本subConductor是不是通孔，如果是的话，那么在生成高斯面时其膨胀的程度要单独处理
    //isVia的默认值是false;
    bool isVia;

    SubConductor();
    SubConductor( int x1, int x2, int y1, int y2, int z1, int z2 );
    SubConductor( const Rectangle &rect );
    SubConductor( const Rectangle &rect,  int fatherCondID );
    SubConductor( const Rectangle &rect,  int fatherCondID , bool viaOrNot );


//@author:Dragon
//judge if two subconductors are actually the same one:
    bool operator == (SubConductor &rhs );

};

typedef std::list<SubConductor> SubConductorList;
//**
//* Conductor
//* - data structure
//*   vec.vec.list.Rectangle
//*   layer.dir.list.Rectangle
class Conductor{
public:
    //**
    //* Dir enum
    //* - corresponds to rectangles with
    //*   -x, +x, -y, +y, -z, +z -dir normals
    enum Dir {LEFT, RIGHT, BACK, FRONT, BOTTOM, TOP};


    static const int nDir = 6;

    int nMetal;
    int nVia;
    int nLayer;

    //it denotes the unique index of the conductor, it is assigned when creating the conductor
    int ID;

    //**
    //* main data
    //* - layer is hierarchized by layer and rect directions.
    //*   This is because
    //*   - It is easier to render
    //*   - It is easier to construct basis functions
    std::vector<std::vector<RectangleList> > layer;

    //@author:Dragon
    //add a subConductorList to express the rectangles included in a conductor , and it also changes when combine two conductors.
   
    SubConductorList subConductorList;

    //@author:Dragon
    //add a gaussianSurfaceList to contain the GaussianSurfaces of all the subConductors.
    GaussianSurfaceList gaussianSurfaceList;
  
    //*
    explicit Conductor();
    explicit Conductor(int nMetal, int nVia);
    explicit Conductor(int nMetal, int nVia, int indexID);

    //**
    //* operator +=
    //* - appended by rhs
    // Conductor& operator+= ( const Conductor& rhs ) throw ( ConductorLayerNotCompatibleError );
    //@author:Dragon
    //because I need to read the subConductorList of conductor rhs, so the type of rhs cannot be const, or there will be error.
    Conductor& operator+= (  Conductor& rhs ) throw ( ConductorLayerNotCompatibleError );

    //**
    //* isContaining
    //* - Check if the 2D +z rect is contained in this conductor regardless of the z-coordinate
    bool isContaining(const Rectangle rect2d, const int layerIndex);

    //* generateVia
    //* - ASSUMPTION: no metal rect has a corner inside the via
    //* - Use SHRINK_VIA macro variable to toggle whether to modify via size or not
    //*   (for the test example generation purpose)
    //* - If isDecomposed is true, then non-overlapping panel decomposition is performed.
    //* - If isDecomposed is false, then overlapping rect is generated for via tops and bottoms.
    void generateVia(
            Rectangle        rect,      // modified only when SHRINK_VIA is defined
            const int        viaIndex,
            const int *const *viaDef,
            const int *const *viaConnect,
            const bool       isDecomposed = false );
    void print() const;

    //**
    //* checkSelftOverlapping
    //* - return true if some rectangles overlap
    //* - only check rects of the same directions within a layer and
    //*   via/metal interfaces
    bool checkSelfOverlapping( const int *const *viaConnect);

    //**
    //* checkZeroAreaRectangele
    //* - return true if this has any rectangle with zero area
    bool checkZeroAreaRectangle();
};


//****
//*
//* Conductor collection
//*

//**
//* ConductorList
//@author:dragon
//add another function into it: get the iterator of the last element, this may cause some errors!!----------This does cause errors.....
// typedef std::list<Conductor> ConductorList;
class ConductorList:public std::list<Conductor>{
public:
    ConductorList::iterator lastIt(){
        ConductorList::iterator iter1,iter2;
        iter1=iter2=this->begin();
        while(++iter2!=this->end()){
            iter1++;
        }
        return iter1;
    }

};

//**
//* LayeredConductorList
typedef std::vector<ConductorList> LayeredConductorList;


//**
//* RectangleGL
//* - Tranform integer-defined Rectangle into single-floating-point
//*   defined RectangleGL for the OpenGL rendering purpose
class RectangleGL
{
public:
    enum ShapeDir  {X_DECAY=0, Y_DECAY=1, Z_DECAY=2, FLAT_SHAPE=3};
    enum ShapeType {FLAT_TYPE, ARCH_TYPE, SIDE_TYPE};

    float   xn;
    float   yn;
    float   zn;

    float   x1;
    float   x2;
    float   y1;
    float   y2;
    float   z1;
    float   z2;

    ShapeType   shapeType;
    ShapeDir    shapeDir;
    float       shapeNormalDistance;
    float       shapeShift;


    RectangleGL();
    RectangleGL(const Rectangle & rect, float unit);

    //**
    //* operation
    bool isOverlappingProjection(const RectangleGL &rect) const;
    bool isOverlapping(const RectangleGL &rect) const;
    bool isOverlappingOrEdgeNeighboring(const RectangleGL &rect) const;
    bool isEmpty() const;
    bool operator== (const RectangleGL &rect) const;
    bool isCoincidental(const RectangleGL &rect, const float margin=0) const;
    bool isContaining(const RectangleGL &rect) const;

    //**
    //* intersectProjection
    //* - intersect the projection of rect onto this
    RectangleGL intersectProjection(const RectangleGL &rect) const;

    //**
    //* intersectArch
    //* - assume same elevation in the normal direction (no check)
    //* - intersect this arch with flat rect
    //* - only consider the case when the edge of decaying head is contained in the flat rect
    RectangleGL intersectArchOnFlat(const RectangleGL &flat) const;



    //**
    //* print
    void print() const;
    void printCapletLine(std::ostream &out) const;
    void printCapletLineFlat(std::ostream &out) const;
};

//**
//* RectangleGLList
class RectangleGLList : public std::list<RectangleGL>
{
public:
    typedef std::list<RectangleGLList::iterator> IteratorList;

    explicit RectangleGLList(
            const allocator_type    &allo = allocator_type());
    explicit RectangleGLList(
            size_type               n,
            const RectangleGL       &value = RectangleGL(),
            const allocator_type    &allo = allocator_type());
    RectangleGLList(
            iterator                first,
            iterator                last,
            const allocator_type    &allo = allocator_type());
    RectangleGLList( const RectangleGLList& rectList);

    //**
    //* mergeProjection
    //* - self merge
    void mergeProjection();
    RectangleGLList mergeProjectionReturn();
    void mergeProjection1_1(const float projectionMergeDistance);

    //**
    //* insertOverlappingRectangleGL
    //* - perform projection of rect onto this list
    //* - find and insert the intersecting RectangleGLs
    //* - return the list of iterators to the inserted RectangleGL
    IteratorList
    insertProjectedOverlappingRectangleGL(const RectangleGL& rect, const float distance);

    void absorbCommonSupport();
    void absorbCommonSupport(IteratorList &itList);

    void removeBadProjection(float margin);

    //* Ver1.0 obsolete
    //void markCommonSupport();


    //****
    //*
    //* DEBUG tools

    //**
    //* hasOverlappingRectngle()
    bool hasCommonSupport() const;
};


//**
//* DirRectangleGLList
//* - Dir: 0 to 5 for rect outer normal dir -x, +x, -y, +y, -z, +z
typedef std::vector<RectangleGLList> DirRectangleGLList;

//**
//* LayeredDirRectangleGLList
//* - Layer index from each metal to each via in order, nLayer == nMetal + nVia
typedef std::vector<DirRectangleGLList> LayeredDirRectangleGLList;

//**
//* ConductorFloat
class ConductorFP
{
public:
    static const unsigned nDir = 6;
    int nMetal;
    int nVia;
    int nLayer;
    LayeredDirRectangleGLList layer;

    explicit ConductorFP();
    explicit ConductorFP(const int nMetal, const int nVia);
    ConductorFP(const Conductor &cond, const float unit);

    unsigned size() const;


};

class ConductorFPList : public std::list<ConductorFP>
{
public:
    explicit ConductorFPList( const allocator_type &allo = allocator_type());
    explicit ConductorFPList(
            size_type n,
            const ConductorFP &value = ConductorFP(),
            const allocator_type &allo = allocator_type());
    ConductorFPList(
            ConductorFPList::iterator first,
            ConductorFPList::iterator last,
            const allocator_type &allo = allocator_type());
    ConductorFPList( const ConductorFPList &condFGList );
    ConductorFPList( const ConductorList &condList, const float unit);

    //**
    //* constructFrom
    //* - clear first
    //* - copy all contents from condList scaled by unit
    void constructFrom( const ConductorList &condList, const float unit );
};


#endif // GDSGEOMETRY_H

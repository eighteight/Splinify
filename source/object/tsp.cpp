#include "c4d.h"
#include "c4d_symbols.h"
#include "c4d_tools.h"
#include "ge_dynamicarray.h"
#include "octsplinify.h"
#include "kd_tree.h"
#include <vector>
#include <string>

typedef std::pair<LONG,Real> longRealPair;
bool comparator ( const longRealPair& l, const longRealPair& r)
{ return l.first < r.first; }

class TSPData : public ObjectData
{
	private:
		LineObject *PrepareSingleSpline(BaseObject *generator, BaseObject *op, Matrix *ml, HierarchyHelp *hh, Bool *dirty);
		void Transform(PointObject *op, const Matrix &m);
		void DoRecursion(BaseObject *op, BaseObject *child, GeDynamicArray<Vector> &points, Matrix ml);
		Random rng;
	public:
		virtual BaseObject* GetVirtualObjects(BaseObject *op, HierarchyHelp *hh);
		virtual Bool Init(GeListNode *node);

		static NodeData *Alloc(void) { return gNew TSPData; }
    
    BaseObject *GetVirtualObjects1(BaseObject *op, HierarchyHelp *hh);
};

void TSPData::Transform(PointObject *op, const Matrix &m)
{
	Vector	*padr=op->GetPointW();
	LONG	pcnt=op->GetPointCount(),i;
	
	for (i=0; i<pcnt; i++)
		padr[i]*=m;
	
	op->Message(MSG_UPDATE);
}

Bool TSPData::Init(GeListNode *node)
{	
	BaseObject		*op   = (BaseObject*)node;
	BaseContainer *data = op->GetDataInstance();

	data->SetReal(CTTSPOBJECT_MAXSEG,3.);
	data->SetBool(CTTSPOBJECT_REL,TRUE);
	return TRUE;
}

void TSPData::DoRecursion(BaseObject *op, BaseObject *child, GeDynamicArray<Vector> &points, Matrix ml)
{
	BaseObject *tp;
	if (child){
		tp = child->GetDeformCache();
		ml = ml * child->GetMl();
		if (tp){
			DoRecursion(op,tp,points,ml);
		}
		else{
			tp = child->GetCache(NULL);
			if (tp){
				DoRecursion(op,tp,points,ml);
			}
			else{

				if (!child->GetBit(BIT_CONTROLOBJECT)){
					if (child->IsInstanceOf(Opoint)){
						PointObject * pChild = ToPoint(child);
						LONG pcnt = pChild->GetPointCount();
						const Vector *childVerts = pChild->GetPointR();
						for(LONG i=0;i<pcnt;i++){
							points.Push(childVerts[i] * ml);
						}
					}
				}
			}
		}
		for (tp = child->GetDown(); tp; tp=tp->GetNext()){
			DoRecursion(op,tp,points,ml);
		}
	}
}

BaseObject *TSPData::GetVirtualObjects1(BaseObject *op, HierarchyHelp *hh)
{
  
	BaseObject *orig = op->GetDown();

	if (!orig) return NULL;

    BaseObject *origSibling = orig->GetNext();
    if (!origSibling) return NULL;


	Bool			dirty=FALSE;
	Matrix			ml, ml1;
	BaseObject *child = op->GetAndCheckHierarchyClone(hh,orig,HIERARCHYCLONEFLAGS_ASPOLY,&dirty,NULL,FALSE);
    BaseObject *sibling = op->GetAndCheckHierarchyClone(hh,origSibling,HIERARCHYCLONEFLAGS_ASPOLY,&dirty,NULL,FALSE);

    if (!sibling) return NULL;

	BaseThread    *bt=hh->GetThread();
	BaseContainer *data = op->GetDataInstance();

	if (!dirty) dirty = op->CheckCache(hh);					
	if (!dirty) dirty = op->IsDirty(DIRTYFLAGS_DATA);		
	if (!dirty) return op->GetCache(hh);

    BaseObject    	*main = BaseObject::Alloc(Onull);
	GeDynamicArray<Vector> childPoints;
	GeDynamicArray<Vector> siblingPoints;
	StatusSetBar(0);
	StatusSetText("Collecting Points");
    //sibling->InsertUnder(child);
	DoRecursion(op,child,childPoints, ml);
	
	StatusSetBar(5);

	rng.Init(1244);
	KDNode *kdTree;
	buildKDTree(childPoints, &kdTree, rng);

    //sibling->Remove();

    PointObject * pChild = ToPoint(child);
    LONG pcnt = pChild->GetPointCount();
    const Vector *childVerts = pChild->GetPointR();
    ml1 = ml1 * child->GetMl();
    std::vector<longRealPair> childPairs(pcnt);
    for(LONG i=0;i<pcnt;i++){
        childPairs[i] = longRealPair(i, 3.0);
        childPoints.Push(childVerts[i] * ml1);
    }
    
	pChild = ToPoint(sibling);
    LONG spcnt = pChild->GetPointCount();
    childVerts = pChild->GetPointR();
    ml = ml * sibling->GetMl();
    std::vector<longRealPair> sibPairs (spcnt);
    for(LONG i=0;i<spcnt;i++){
        siblingPoints.Push(childVerts[i] * ml);
    }
    
    for(LONG i=0;i<spcnt;i++){
        Real dist = MAXREALr;
        LONG minInd;
        for (int j = 0; j<childPoints.GetCount(); j++){
            Real newDist = (siblingPoints[i] - childPoints[j]).GetLength();
            if (newDist<dist){
                for (int k=0; i < i; k++){
                    if (sibPairs[k].first == j && sibPairs[k].second>dist) {
                        dist = newDist;
                        minInd = j;
                    }
                }
            }
        }
        sibPairs[i] = longRealPair(minInd, dist);
    }
    
    sort(sibPairs.begin(), sibPairs.end(), [](longRealPair const& a, longRealPair const& b)
    {
        if (a.second < b.second) return true;
        if (a.second > b.second) return false;
        
        return false;
    });

	if(spcnt > 0){
		GeDynamicArray<LONG> pointList(pcnt);
		GeDynamicArray<LONG> path(pcnt);
		LONG currentPoint = 0;
		for(LONG i=0;i<pcnt;i++){
			pointList[i] = 1;
		}
		path[0] = 0;
		pointList[0] = 0;
		StatusSetText("Connecting Points");
		Real dist;
		for(LONG i=0;i<spcnt;i++){
			dist = MAXREALr;
            LONG minInd;
            Vector p1 = siblingPoints[i];
            for (int j = 0; j<childPoints.GetCount(); j++){
                Real newDist = (p1 - childPoints[j]).GetLength();
                if (newDist<dist){
                    dist = newDist;
                    minInd = j;
                }
            }

            Vector *padr;
            SplineObject	*spline=SplineObject::Alloc(2,SPLINETYPE_LINEAR); //http://www.microbion.co.uk/graphics/c4d/create_plugins4b.htm
            if (!spline) continue;
            spline->GetDataInstance()->SetBool(SPLINEOBJECT_CLOSED, FALSE);
            padr = spline->GetPointW();

            Vector p2 = childPoints[minInd];
            padr[0] = p1;
            padr[1] = p2;
			spline->Message(MSG_UPDATE);
            spline->SetName(child->GetName());
            spline->InsertUnder(main);
            if(i % 20 == 0){
				StatusSetBar(10 + (90*i)/pcnt);
				if (bt && bt->TestBreak()){
					pcnt = i;
					break;
				}
			}
        }
	}
	GeFree(kdTree);
	StatusClear();

	return main;
}

BaseObject *TSPData::GetVirtualObjects(BaseObject *op, HierarchyHelp *hh)
{
    
	BaseObject *orig = op->GetDown();
    
	if (!orig) return NULL;
    
    BaseObject *origSibling = orig->GetNext();
    if (!origSibling) return NULL;
    
    // start new list
	op->NewDependenceList();

	Bool			dirty=FALSE;
	Matrix			ml;
	BaseObject *child = op->GetAndCheckHierarchyClone(hh,orig,HIERARCHYCLONEFLAGS_ASPOLY,&dirty,NULL,FALSE);
    BaseObject *sibling = op->GetAndCheckHierarchyClone(hh,origSibling,HIERARCHYCLONEFLAGS_ASPOLY,&dirty,NULL,FALSE);
    
	BaseThread    *bt=hh->GetThread();
	BaseContainer *data = op->GetDataInstance();
	Real maxSeg = data->GetReal(CTTSPOBJECT_MAXSEG,3.);
	Bool relativeMaxSeg  = data->GetBool(CTTSPOBJECT_REL,TRUE);
    
	if (!dirty) dirty = op->CheckCache(hh);
	if (!dirty) dirty = op->IsDirty(DIRTYFLAGS_DATA);
    // mark child objects as processed
	op->TouchDependenceList();
	if (!dirty) return op->GetCache(hh);
    
    BaseObject    	*main = BaseObject::Alloc(Onull);
	GeDynamicArray<Vector> childPoints;
	GeDynamicArray<Vector> siblingPoints;
	StatusSetBar(0);
	StatusSetText("Collecting Points");
	DoRecursion(op,child,childPoints, ml);
	DoRecursion(op,sibling,siblingPoints, ml);
    
	StatusSetBar(5);
    
	rng.Init(1244);
	KDNode *kdTree;
	buildKDTree(childPoints, &kdTree, rng);

	LONG pcnt = siblingPoints.GetCount();
	if(pcnt > 0){
		GeDynamicArray<LONG> pointList(pcnt);

		LONG currentPoint = 0;
		for(LONG i=0;i<pcnt;i++){
			pointList[i] = 1;
		}

		pointList[0] = 0;

		StatusSetText("Connecting Points");
		Real dist;
		for(LONG i=0;i<pcnt;i++){
			dist = -1.;
			LONG closestPoint = kdTree->getNearestNeighbor(childPoints,siblingPoints[i],pointList, dist, 0);
			if(closestPoint == -1){
				GePrint("error finding neighbor");
				pcnt = i-1;
				break;
			}
            
			pointList[closestPoint] = 0;
            
            Vector *padr;
            SplineObject	*spline=SplineObject::Alloc(2,SPLINETYPE_LINEAR); //http://www.microbion.co.uk/graphics/c4d/create_plugins4b.htm
            if (!spline) continue;
            spline->GetDataInstance()->SetBool(SPLINEOBJECT_CLOSED, FALSE);
            padr = spline->GetPointW();
            
            Vector p2 = childPoints[closestPoint];
            padr[0] = siblingPoints[i];
            padr[1] = p2;
			spline->Message(MSG_UPDATE);
            spline->SetName(child->GetName());
            spline->InsertUnder(main);
            if(i % 20 == 0){
				StatusSetBar(10 + (90*i)/pcnt);
				if (bt && bt->TestBreak()){
					pcnt = i;
					break;
				}
			}
		}
	}
	GeFree(kdTree);
	StatusClear();
	return main;
Error:
	BaseObject::Free(child);
    BaseObject::Free(child);
	return NULL;
}

// be sure to use a unique ID obtained from www.plugincafe.com
#define ID_TSPOBJECT 1030923

Bool RegisterTSP(void)
{
	return RegisterObjectPlugin(ID_TSPOBJECT,GeLoadString(IDS_TSP),OBJECT_GENERATOR|OBJECT_INPUT,TSPData::Alloc,"Octsplinify",AutoBitmap("tsp.tif"),0);
}

#include "c4d.h"
#include "c4d_symbols.h"
#include "c4d_tools.h"
#include "ge_dynamicarray.h"
#include "octsplinify.h"
#include "kd_tree.h"
#include <vector>
#include <string>

#include <iostream>

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
        Bool isCalculated;
	public:
		virtual BaseObject* GetVirtualObjects(BaseObject *op, HierarchyHelp *hh);
		virtual Bool Init(GeListNode *node);

		static NodeData *Alloc(void) { return gNew TSPData; }
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

	data->SetReal(CTTSPOBJECT_MAXSEG,30.);
	data->SetBool(CTTSPOBJECT_REL,TRUE);
    isCalculated = FALSE;
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

BaseObject *TSPData::GetVirtualObjects(BaseObject *op, HierarchyHelp *hh)
{
    
	BaseObject *orig = op->GetDown();
    
	if (!orig) return NULL;
    
    // start new list
	op->NewDependenceList();
    
	// check cache for validity and check master object for changes
	Bool dirty = op->CheckCache(hh) || op->IsDirty(DIRTYFLAGS_DATA);
    BaseObject* pp=NULL;
	// for each child
    0;

    GeDynamicArray<BaseObject*> children;
	for (pp=orig; pp; pp=pp->GetNext()) {
        children.Push(op->GetHierarchyClone(hh,pp,HIERARCHYCLONEFLAGS_ASPOLY,FALSE,NULL));
	}
    LONG child_cnt = children.GetCount();
    
	// no child objects found
	if (child_cnt == 0) return NULL;
    
	// if child list has been modified
	if (!dirty) dirty = !op->CompareDependenceList();
    
	// mark child objects as processed
	op->TouchDependenceList();
    
	// if no change has been detected, return original cache
	if (!dirty) return op->GetCache(hh);
    std::cout<<"no cache"<<std::endl;

    BaseContainer *data = op->GetDataInstance();
    Real maxSeg = data->GetReal(CTTSPOBJECT_MAXSEG,30.);
	Bool relativeMaxSeg  = data->GetBool(CTTSPOBJECT_REL,TRUE);
    BaseThread    *bt=hh->GetThread();
    BaseObject* main = BaseObject::Alloc(Onull);
    isCalculated = TRUE;
	StatusSetBar(0);
    StatusSetText("Collecting Points");
    GeDynamicArray<KDNode*> trees(child_cnt);
    GeDynamicArray<GeDynamicArray<Vector> > chldPoints(child_cnt);
    GeDynamicArray<GeDynamicArray<LONG> > validPoints(child_cnt);
    rng.Init(1244);
    
    LONG maxPointCnt = 0;
    for (int k=0; k < child_cnt; k++){
        Matrix ml;
        DoRecursion(op,children[k],chldPoints[k], ml);
        KDNode *kdTree;
        buildKDTree(chldPoints[k], &kdTree, rng);
        trees[k] = kdTree;
        validPoints[k] = GeDynamicArray<LONG>(chldPoints[k].GetCount());
        for (LONG i=0; i < chldPoints[k].GetCount(); i++){
            validPoints[k][i] = 1;
        }
        if (maxPointCnt< chldPoints[k].GetCount()){
            maxPointCnt = chldPoints[k].GetCount();
        }
    }

    StatusSetBar(5);
    StatusSetText("Connecting Points");

    for (LONG i = 0; i < maxPointCnt; i++){
        GeDynamicArray<Vector> splinePoints;
        LONG goodCnt = 0;
        for (int k=0; k < child_cnt-1; k++){
            GePrint(children[k]->GetName()+" ->"+children[k+1]->GetName());
            if (chldPoints[k].GetCount() < i || chldPoints[k+1].GetCount()<i) continue;

            validPoints[k][0] = 0;
            
            Vector nxt = chldPoints[k+1][i];

            Real dist = -1.;
            LONG closestPointIndx = trees[k]->getNearestNeighbor(chldPoints[k], nxt, validPoints[k], dist, 0);
            if(closestPointIndx == -1){
                GePrint("error finding neighbor");
                continue;
            }

            std::cout<<"d: "<<dist<<std::endl;
            if (dist> maxSeg || dist < 0.01) {
                continue;
            }
            goodCnt++;
            validPoints[k][closestPointIndx] = 0;
            
            Vector closest = chldPoints[k][closestPointIndx];
            splinePoints.Push(closest);
            splinePoints.Push(nxt);
        }
        
        if (splinePoints.GetCount() == 0) continue;
        
        SplineObject	*spline=SplineObject::Alloc(splinePoints.GetCount(),SPLINETYPE_LINEAR);
        if (!spline) continue;
        spline->GetDataInstance()->SetBool(SPLINEOBJECT_CLOSED, FALSE);
        Vector *padr = spline->GetPointW();
        for (LONG l=0;l<splinePoints.GetCount();l++){
            GePrint(RealToString(splinePoints[l].x)+" "+RealToString(splinePoints[l].y)+" "+RealToString(splinePoints[l].z));
            padr[l] = splinePoints[l];
        }
        spline->InsertUnder(main);
        spline->Message(MSG_UPDATE);
        
        if(i % 20 == 0){
            StatusSetBar(10 + (90*i)/maxPointCnt);
            if (bt && bt->TestBreak()){
                //pcnt = i;
                break;
            }
        }
        GePrint(LongToString(goodCnt));
    }
    
    for (int k=0; k<child_cnt; k++){
        GeFree(trees[k]);
    }
    
    main->Message(MSG_UPDATE);
	StatusClear();
	return main;
Error:
//	BaseObject::Free(childs[0]);
//    BaseObject::Free(childs[1]);
	return NULL;
}

// unique ID obtained from www.plugincafe.com
#define ID_TSPOBJECT 1030923

Bool RegisterTSP(void)
{
	return RegisterObjectPlugin(ID_TSPOBJECT,GeLoadString(IDS_TSP),OBJECT_GENERATOR|OBJECT_INPUT,TSPData::Alloc,"Octsplinify",AutoBitmap("tsp.tif"),0);
}

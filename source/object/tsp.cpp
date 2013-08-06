#include "c4d.h"
#include "c4d_symbols.h"
#include "c4d_tools.h"
#include "ge_dynamicarray.h"
#include "octsplinify.h"
#include "kd_tree.h"
#include <vector>
#include <string>

#include <iostream>


#define MAXTARGETS 5
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
    LONG child_cnt =0;

    BaseObject* children[MAXTARGETS];
	for (pp=orig; pp; pp=pp->GetNext()) {
        children[child_cnt++]=op->GetHierarchyClone(hh,pp,HIERARCHYCLONEFLAGS_ASPOLY,FALSE,NULL);
        if ((child_cnt)== MAXTARGETS) break;
	}
    
	// no child objects found
	if (!child_cnt) return NULL;
    
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
    rng.Init(1244);
    for (int k=0; k < child_cnt; k++){
        Matrix			ml;
        GeDynamicArray<Vector> childPoints;
        DoRecursion(op,children[k],chldPoints[k], ml);
        KDNode *kdTree;
        buildKDTree(chldPoints[k], &kdTree, rng);
        trees[k] = kdTree;
    }

    GeDynamicArray<GeDynamicArray<LONG> > valdPoints(child_cnt);
    
    for(LONG k=0;k<child_cnt;k++){
        valdPoints[k] = GeDynamicArray<LONG>(chldPoints[k].GetCount());
        for (LONG i=0; i < chldPoints[k].GetCount(); i++){
            valdPoints[k][i] = 1;
        }
    }
    
    StatusSetBar(5);
    StatusSetText("Connecting Points");
    for (int k=0; k < 1; k++){
        LONG pcnt = chldPoints[k+1].GetCount();
        LONG goodCnt = 0;
        if(pcnt > 0){
            GeDynamicArray<LONG> validPoints(pcnt);

            for(LONG i=0;i<pcnt;i++){
                validPoints[i] = 1;
            }

            validPoints[0] = 0;

            Real dist;
            for(LONG i=0; i < pcnt; i++){
                dist = -1.;
                LONG closestPointIndx = trees[k]->getNearestNeighbor(chldPoints[k],chldPoints[k+1][i],validPoints, dist, 0);
                if(closestPointIndx == -1){
                    GePrint("error finding neighbor");
                    continue;
                }
                
                std::cout<<"d: "<<dist<<std::endl;
                if (dist> maxSeg || dist < 0.01) {
                    continue;
                }
                goodCnt++;
                validPoints[closestPointIndx] = 0;
                
                Vector *padr;
                SplineObject	*spline=SplineObject::Alloc(2,SPLINETYPE_LINEAR);
                if (!spline) continue;
                spline->GetDataInstance()->SetBool(SPLINEOBJECT_CLOSED, FALSE);
                padr = spline->GetPointW();
                
                Vector p2 = chldPoints[k][closestPointIndx];
                padr[0] = chldPoints[k+1][i];
                padr[1] = p2;
                spline->SetName(children[0]->GetName());
                spline->InsertUnder(main);
                main->Message(MSG_UPDATE);
                if(i % 20 == 0){
                    StatusSetBar(10 + (90*i)/pcnt);
                    if (bt && bt->TestBreak()){
                        pcnt = i;
                        break;
                    }
                }
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

// be sure to use a unique ID obtained from www.plugincafe.com
#define ID_TSPOBJECT 1030923

Bool RegisterTSP(void)
{
	return RegisterObjectPlugin(ID_TSPOBJECT,GeLoadString(IDS_TSP),OBJECT_GENERATOR|OBJECT_INPUT,TSPData::Alloc,"Octsplinify",AutoBitmap("tsp.tif"),0);
}

#include "c4d.h"
#include "c4d_symbols.h"
#include "c4d_tools.h"
#include "customgui_inexclude.h"
#include "ge_dynamicarray.h"
#include "octsplinify.h"
#include "kd_tree.h"
#include <vector>
#include <string>

#include <iostream>

typedef std::pair<LONG,Real> longRealPair;
bool comparator ( const longRealPair& l, const longRealPair& r)
{ return l.first < r.first; }

class SplinifyData : public ObjectData
{
	private:
		LineObject *PrepareSingleSpline(BaseObject *generator, BaseObject *op, Matrix *ml, HierarchyHelp *hh, Bool *dirty);
		void Transform(PointObject *op, const Matrix &m);
		void DoRecursion(BaseObject *op, BaseObject *child, GeDynamicArray<Vector> &points, Matrix ml);
        Random rng;
        Bool isCalculated;
        GeDynamicArray<GeDynamicArray<Vector> > splineAtPoint;
        LONG prvsFrame = 0;
	public:
		virtual BaseObject* GetVirtualObjects(BaseObject *op, HierarchyHelp *hh);
		virtual Bool Init(GeListNode *node);

		static NodeData *Alloc(void) { return gNew SplinifyData; }
};

void SplinifyData::Transform(PointObject *op, const Matrix &m)
{
	Vector	*padr=op->GetPointW();
	LONG	pcnt=op->GetPointCount(),i;
	
	for (i=0; i<pcnt; i++)
		padr[i]*=m;
	
	op->Message(MSG_UPDATE);
}

Bool SplinifyData::Init(GeListNode *node)
{	
	BaseObject		*op   = (BaseObject*)node;
	BaseContainer *data = op->GetDataInstance();

	data->SetReal(CTTSPOBJECT_MAXSEG,30.);
	data->SetBool(CTTSPOBJECT_REL,TRUE);
    data->SetLong(SPLINEOBJECT_INTERPOLATION,SPLINEOBJECT_INTERPOLATION_ADAPTIVE);
    isCalculated = FALSE;
    GePrint("Splinify by http://twitter.com/eight_io for Cinema 4D r14");
	return TRUE;
}

void SplinifyData::DoRecursion(BaseObject *op, BaseObject *child, GeDynamicArray<Vector> &points, Matrix ml)
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

BaseObject *SplinifyData::GetVirtualObjects(BaseObject *op, HierarchyHelp *hh)
{
    BaseDocument *doc = op->GetDocument();
    LONG crntFrame = doc->GetTime().GetFrame(doc->GetFps());
    
    BaseContainer *data = op->GetDataInstance();
    LONG delta = data->GetLong(CTTSPOBJECT_WINDOW,1);
    LONG strtFrame = crntFrame - delta;
    strtFrame = strtFrame<0?0:strtFrame;
    // start new list
	op->NewDependenceList();
    
	// check cache for validity and check master object for changes
	Bool dirty = op->CheckCache(hh) || op->IsDirty(DIRTYFLAGS_DATA);
    

    InExcludeData* traceElements = (InExcludeData*)data->GetCustomDataType(TRACE_ELEMENTS, CUSTOMDATATYPE_INEXCLUDE_LIST);
    GeDynamicArray<BaseObject*> children;
    if(traceElements && traceElements->GetObjectCount()>0) {
        for(int i=0; i<traceElements->GetObjectCount(); ++i) {
            BaseObject* pp = (BaseObject*)traceElements->ObjectFromIndex(op->GetDocument(),i);
            if (pp) {
                BaseObject* chld = NULL;
                LONG trck = 0;
                for (chld=pp->GetDownLast(); chld; chld=chld->GetPred()) {
                    if (trck > strtFrame && trck<= crntFrame){
                        children.Push(op->GetHierarchyClone(hh,chld,HIERARCHYCLONEFLAGS_ASPOLY,FALSE,NULL));
                    }
                    trck++;
                }
            }
            if (i == 0) break; //consider only one child
        }
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
    
    GePrint("NO CACHE");

    Real maxSeg = data->GetReal(CTTSPOBJECT_MAXSEG,30.);
	Bool relativeMaxSeg  = data->GetBool(CTTSPOBJECT_REL,TRUE);

    LONG splineInterpolation = data->GetLong(SPLINEOBJECT_INTERPOLATION);

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
        if (maxPointCnt < chldPoints[k].GetCount()){
            maxPointCnt = chldPoints[k].GetCount();
        }
    }
    
    for (int k=0; k < child_cnt; k++){
        validPoints[k] = GeDynamicArray<LONG>(maxPointCnt);
        validPoints[k].Fill(0,maxPointCnt,1);
    }

    StatusSetBar(5);
    StatusSetText("Connecting Points");

    Real distMin = MAXREALr;
    Real distMax = 0.;

    splineAtPoint.FreeArray();
    
    LONG shift = crntFrame - prvsFrame;
    shift = shift < 0? 0: shift;
    
    for (LONG i = 0; i < splineAtPoint.GetCount(); i++){
        splineAtPoint[i].Shift(0, -shift);
        splineAtPoint[i].ReSize(splineAtPoint[i].GetCount() - shift);
    }

    for (LONG i = 0; i < maxPointCnt; i++){
        if (splineAtPoint.GetCount() == i) {
            splineAtPoint.Push(GeDynamicArray<Vector>());
        }
        
        for (int k=0; k < child_cnt-1; k++){
            if (chldPoints[k].GetCount() < i || chldPoints[k+1].GetCount()<i || validPoints[k].GetCount()<i) {
                continue;
            }

            validPoints[k][0] = 0;
            Vector queryPoint = splineAtPoint[i].GetCount() == 0? chldPoints[k][i]:splineAtPoint[i][splineAtPoint[i].GetCount()-1];

            Real dist = -1.;
            LONG closestIndx = trees[k+1]->getNearestNeighbor(chldPoints[k+1], queryPoint, validPoints[k], dist, 0);
            if(closestIndx == -1){
                GePrint("error finding neighbor");
                continue;
            }
            
            distMin = distMin<dist?distMin:dist;
            distMax = distMax>dist?distMax:dist;
            
            if (dist > maxSeg || dist < 0.01) {
                continue;
            }
            validPoints[k][closestIndx] = 0;

            splineAtPoint[i].Push(chldPoints[k+1][closestIndx]);
        }
        
        if (splineAtPoint[i].GetCount() == 0||splineAtPoint[i].GetCount()<10) continue;
        
        SplineObject	*spline=SplineObject::Alloc(splineAtPoint[i].GetCount(),SPLINETYPE_AKIMA);
        if (!spline) continue;
        spline->GetDataInstance()->SetBool(SPLINEOBJECT_CLOSED, FALSE);
        Vector *padr = spline->GetPointW();
        for (LONG l=0;l<splineAtPoint[i].GetCount();l++){
            padr[l] = splineAtPoint[i][l];
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
    }
    GePrint(LongToString(strtFrame)+"-"+LongToString(crntFrame)+ ":"+LongToString(splineAtPoint.GetCount())+" dist: "+RealToString(distMin)+":"+RealToString(distMax));
    
    for (int k=0; k<child_cnt; k++){
        GeFree(trees[k]);
        if (children[k]){
            BaseObject::Free(children[k]);
        }
    }
    main->Message(MSG_UPDATE);
    prvsFrame = crntFrame;
	StatusClear();
	return main;
Error:
    for (int i = 0; i < children.GetCount(); i++){
        BaseObject::Free(children[i]);
    }
	return NULL;
}

// unique ID obtained from www.plugincafe.com
#define ID_SPLINIFYOBJECT 1030923

Bool RegisterTSP(void)
{
	return RegisterObjectPlugin(ID_SPLINIFYOBJECT,GeLoadString(IDS_SPLINIFY),OBJECT_GENERATOR|OBJECT_INPUT,SplinifyData::Alloc,"Octsplinify",AutoBitmap("tsp.tif"),0);
}

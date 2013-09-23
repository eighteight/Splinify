#include "c4d.h"
#include "c4d_symbols.h"
#include "c4d_tools.h"
#include "lib_splinehelp.h"
#include "ge_dynamicarray.h"
#include "octsplinify.h"
#include "kd_tree.h"
#include <vector>
#include <string>

#include <iostream>

// unique ID obtained from www.plugincafe.com
#define ID_SPLINIFYOBJECT 1030923

typedef std::pair<SplineObject*,Real> SplinePair;
bool comparator ( const SplinePair& l, const SplinePair& r){
    return l.second > r.second;
}

class SplinifyData : public ObjectData
{
	private:
		void DoRecursion(BaseObject *op, BaseObject *child, GeDynamicArray<Vector> &points, Matrix ml);
        Random rng;
        Real maxSeg, minSeg;
        Matrix parentMatrix;
        SplineObject* ComputeSpline(BaseThread* bt, GeDynamicArray<GeDynamicArray<Vector> > &objectPoints, GeDynamicArray<KDNode*> &trees, LONG maxPoints, LONG longestPercent, GeDynamicArray<GeDynamicArray<Vector> > &splineAtPoint);
    
	public:
		virtual SplineObject* GetContour(BaseObject *op, BaseDocument *doc, Real lod, BaseThread *bt);
		virtual Bool Init(GeListNode *node);
		static NodeData *Alloc(void) { return gNew SplinifyData; }
};

Bool SplinifyData::Init(GeListNode *node)
{	
	BaseObject		*op   = (BaseObject*)node;
	BaseContainer *data = op->GetDataInstance();

	data->SetReal(CTTSPOBJECT_MAXSEG,1000.);
	data->SetReal(CTTSPOBJECT_MINSEG,0.1);
    data->SetLong(SPLINEOBJECT_INTERPOLATION,SPLINEOBJECT_INTERPOLATION_ADAPTIVE);
    GePrint("Splinify by http://twitter.com/eight_io for Cinema 4D r14");
    
    LONG i;
    GeData param;
    const DescID id = DescID(ID_SPLINIFYOBJECT);
    
    // set up the initial parameter for the Execute function
    i = 0;
    param.SetLong(i);
    op->SetParameter(id, param, DESCFLAGS_SET_0);
    return TRUE;
}

void SplinifyData::DoRecursion(BaseObject *op, BaseObject *child, GeDynamicArray<Vector> &points, Matrix ml) {
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
						for(LONG i=0; i < pcnt; i++){
							points.Push(childVerts[i] * ml * parentMatrix);
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

SplineObject* SplinifyData::GetContour(BaseObject *op, BaseDocument *doc, Real lod, BaseThread *bt){
    BaseContainer *data = op->GetDataInstance();
    BaseObject* parent=(BaseObject*)data->GetLink(CTT_OBJECT_LINK,doc,Obase);
    if (!parent) return NULL;
    
    LONG startObject = data->GetLong(START_FRAME);
    LONG endObject = data->GetLong(END_FRAME);
    
    if (startObject >=endObject) return NULL;
    
    maxSeg = data->GetReal(CTTSPOBJECT_MAXSEG,30.);
    minSeg = data->GetReal(CTTSPOBJECT_MINSEG);

    LONG delta = data->GetLong(OBJECT_SKIP,1);
    delta = delta < 1? 1 : delta;

    GeDynamicArray<BaseObject*> children;
    GeDynamicArray<GeDynamicArray<Vector> > splineAtPoint;
    
    BaseObject* chld = NULL;
    LONG trck = 0;
    for (chld=parent->GetDownLast(); chld; chld=chld->GetPred()) {
        if (trck >= startObject && trck<= endObject && trck % delta == 0){
            children.Push((BaseObject*)chld->GetClone(COPYFLAGS_NO_HIERARCHY|COPYFLAGS_NO_ANIMATION|COPYFLAGS_NO_BITS,NULL));
        }
        trck++;
    }
    
    LONG child_cnt = children.GetCount();
    if (child_cnt < 2) {
        return NULL;
    }

    LONG splineInterpolation = data->GetLong(SPLINEOBJECT_INTERPOLATION);
    LONG longestPercent = data->GetLong(TAKE_LONGEST, 1);
    longestPercent = longestPercent > 100 ? 100: longestPercent;
    
	StatusSetBar(0);
    StatusSetText("Collecting Points");
    GeDynamicArray<KDNode*> trees(children.GetCount());
    GeDynamicArray<GeDynamicArray<Vector> > objectPoints(child_cnt);
    
    rng.Init(1244);

    LONG maxPointCnt = 0;
    parentMatrix = parent->GetMl();
    
    for (int k= 0; k < children.GetCount(); k++){
        GePrint(children[k]->GetName());
        Matrix ml;
        DoRecursion(op,children[k],objectPoints[k], ml);
        KDNode *kdTree;
        buildKDTree(objectPoints[k], &kdTree, rng);
        trees[k] = kdTree;
        if (maxPointCnt < objectPoints[k].GetCount()){
            maxPointCnt = objectPoints[k].GetCount();
        }
    }
    
    SplineObject* parentSpline = ComputeSpline(bt, objectPoints, trees, maxPointCnt, longestPercent,splineAtPoint);

    ModelingCommandData mcd;
    mcd.doc = doc;
    mcd.op = parentSpline;
    
    if(!SendModelingCommand(MCOMMAND_JOIN, mcd)){
        return NULL;
    }
    
    SplineObject* ret = ToSpline(mcd.result->GetIndex(0L));
    
    ret->GetDataInstance()->SetLong(SPLINEOBJECT_INTERPOLATION, splineInterpolation);
    
    for (int k=0; k<child_cnt; k++){
        GeFree(trees[k]);
    }
    
    for (int k=0; k<children.GetCount(); k++){
        if (children[k]){
            BaseObject::Free(children[k]);
        }
    }

    return ret;
Error:
    for (int i = 0; i < children.GetCount(); i++){
        BaseObject::Free(children[i]);
    }
    return NULL;
}

SplineObject* SplinifyData::ComputeSpline(BaseThread* bt, GeDynamicArray<GeDynamicArray<Vector> > &objectPoints, GeDynamicArray<KDNode*> &trees, LONG maxPoints, LONG longestPercent, GeDynamicArray<GeDynamicArray<Vector> > &splinesAtPoint){

    StatusSetBar(5);
    StatusSetText("Connecting Points");
    
    SplineObject* parentSpline = SplineObject::Alloc(0, SPLINETYPE_BSPLINE);
    
    std::vector<SplinePair >splinePairs;

    Real avSplineSize = 0.0, avSplineLength = 0.0;

    if (splinesAtPoint.GetCount() == 0){
        Random r;
        r.Init(43432);
        GeDynamicArray<LONG> validPoints(objectPoints[0].GetCount());
        GeDynamicArray<LONG> usedIndxes(objectPoints[0].GetCount());

        validPoints.Fill(0,objectPoints[0].GetCount(),1);

        for (LONG i = 0; i < objectPoints[0].GetCount(); i++) {

            Vector queryPoint = objectPoints[0][i];
            
            Real dist = -1.;
            LONG closestIndx = trees[1]->getNearestNeighbor(objectPoints[1], queryPoint, validPoints, dist, 0);
            if (closestIndx == -1){
                GePrint("error finding neighbor");
                continue;
            }
            
            if (dist > maxSeg || dist < minSeg) {
                continue;
            }

            if (dist > 0.0){
                validPoints[closestIndx] = 0;
                GeDynamicArray<Vector> rawSpline;
                rawSpline.Push(queryPoint);
                splinesAtPoint.Push(rawSpline);
            }
        }
    }

    GeDynamicArray<GeDynamicArray<LONG> > validPoints(objectPoints.GetCount());
    for (LONG k=0; k < objectPoints.GetCount(); k++){
        validPoints[k] = GeDynamicArray<LONG>(maxPoints);
        validPoints[k].Fill(0,maxPoints,1);
    }

    Real distMin = MAXREALr;
    Real distMax = 0.;
    AutoAlloc<SplineHelp> splineHelp;
    for (LONG i = 0; i < splinesAtPoint.GetCount(); i++){//iterate points
        GeDynamicArray<Vector> queryVector = splinesAtPoint[i];

        for (LONG o=0; o < objectPoints.GetCount()-1; o++){ //iterate objects

            Vector queryPoint = queryVector[queryVector.GetCount()-1];
            
            Real dist = -1.;
            LONG closestIndx = trees[o+1]->getNearestNeighbor(objectPoints[o+1], queryPoint, validPoints[o], dist, 0); //query next object
            if(closestIndx == -1){
                GePrint("error finding neighbor");
                continue;
            }
            
            distMin = distMin < dist ? distMin : dist;
            distMax = distMax > dist ? distMax : dist;
            
            if (dist > maxSeg || dist < 10.0) {
                continue;
            }
            validPoints[o][closestIndx] = 0;
            Vector clsst = objectPoints[o+1][closestIndx];
            if (splinesAtPoint[i].Find(clsst) == NOTOK){
                splinesAtPoint[i].Push(clsst);
            }
        }

        if (splinesAtPoint[i].GetCount() == 0) continue;
        
        SplineObject	*spline=SplineObject::Alloc(splinesAtPoint[i].GetCount(),SPLINETYPE_BSPLINE);
        if (!spline) continue;
        
        spline->GetDataInstance()->SetBool(SPLINEOBJECT_CLOSED, FALSE);
        
        Vector *padr = spline->GetPointW();
        for (LONG l=0;l<splinesAtPoint[i].GetCount();l++){
            padr[l] = splinesAtPoint[i][l];
        }

        splineHelp->InitSpline(spline);
        Real splnLength = splineHelp->GetSplineLength();
        if (splnLength>0.0){
            splinePairs.push_back(SplinePair(spline, splnLength));
            avSplineLength += splnLength;
            avSplineSize += splinesAtPoint[i].GetCount();
        } else {
            SplineObject::Free(spline);
        }

        if(i % 5 == 0){
            StatusSetBar(10 + (90*i)/splinesAtPoint.GetCount());
            if (bt && bt->TestBreak()){
                break;
            }
        }
    }
    
    LONG splnSize = splinePairs.size();
    if (splnSize > 0) {
        LONG limit =  splnSize * longestPercent / 100;
        limit = limit == 0 ? 1 : limit;
        
        std::sort(splinePairs.begin(), splinePairs.end(),comparator);
        
        for (int s = 0; s < limit; s++){
            avSplineLength += splinePairs[s].second;
            splinePairs[s].first->InsertUnder(parentSpline);
        }

        String sizeAvg = splinesAtPoint.GetCount() == 0? "Nan":RealToString(avSplineSize/splinesAtPoint.GetCount());
        
        GePrint("pnts="+sizeAvg+" d="+RealToString(distMin)+":"+RealToString(distMax)+" sl="+RealToString(avSplineLength/avSplineSize));
    }
	StatusClear();
    
    if (splnSize == 0) return NULL;
    return parentSpline;
}


Bool RegisterSplinify(void)
{
	return RegisterObjectPlugin(ID_SPLINIFYOBJECT,GeLoadString(IDS_SPLINIFY),OBJECT_GENERATOR|OBJECT_INPUT|OBJECT_ISSPLINE,SplinifyData::Alloc,"Octsplinify",AutoBitmap("tsp.tif"),0);
}

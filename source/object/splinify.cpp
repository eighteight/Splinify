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
    return l.first < r.first;
}

class SplinifyData : public ObjectData
{
	private:
		void Transform(PointObject *op, const Matrix &m);
		void DoRecursion(BaseObject *op, BaseObject *child, GeDynamicArray<Vector> &points, Matrix ml);
        Random rng;
        GeDynamicArray<GeDynamicArray<Vector> > splineAtPoint;
        LONG prvsFrame = 0, oldFrame;
        SplineObject* GetSpline(BaseObject* op, BaseThread* bt, BaseDocument *doc, GeDynamicArray<BaseObject*> &children, Real maxSeg, LONG delta, LONG splinePercentage);
    
	public:
		virtual SplineObject* GetContour(BaseObject *op, BaseDocument *doc, Real lod, BaseThread *bt);
		virtual Bool Init(GeListNode *node);
		static NodeData *Alloc(void) { return gNew SplinifyData; }
    
        virtual Bool AddToExecution(BaseObject *op, PriorityList *list);
        virtual EXECUTIONRESULT Execute(BaseObject *op, BaseDocument *doc, BaseThread *bt, LONG priority, EXECUTIONFLAGS flags);
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
    GePrint("Splinify by http://twitter.com/eight_io for Cinema 4D r14");
    
    LONG i;
    GeData param;
    const DescID id = DescID(ID_SPLINIFYOBJECT);
    
    // 'oldframe' is a class-level variable that will let us check if the current frame is different
    // from the one when the Execute() function was last called
    oldFrame = 0;
    
    // set up the initial parameter for the Execute function
    i = 0;
    param.SetLong(i);
    op->SetParameter(id, param, DESCFLAGS_SET_0);
    return TRUE;
}

Bool SplinifyData::AddToExecution(BaseObject *op, PriorityList *list)
{
    list->Add(op, EXECUTIONPRIORITY_GENERATOR, EXECUTIONFLAGS_ANIMATION);
    return TRUE;
}

EXECUTIONRESULT SplinifyData::Execute(BaseObject *op, BaseDocument *doc, BaseThread *bt, LONG priority, EXECUTIONFLAGS flags)
{
    GeData param;
    const DescID id = DescID(ID_SPLINIFYOBJECT);
    LONG i, fps, frame;
    
    fps = doc->GetFps();
    frame = doc->GetTime().GetFrame(fps);
    
    // ensure that this only happens once per frame or GetContour may be called multiple times
    if(frame == 0 || frame != oldFrame)
    {
        // hack to force GetContour to update
        op->GetParameter(id, param, DESCFLAGS_GET_0);
        i = param.GetLong();
        i++;
        param.SetLong(i);
        op->SetParameter(id, param, DESCFLAGS_SET_0);
        oldFrame = frame;
    }
    
    return EXECUTIONRESULT_OK;
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

SplineObject* SplinifyData::GetContour(BaseObject *op, BaseDocument *doc, Real lod, BaseThread *bt){
    BaseContainer *data = op->GetDataInstance();
    BaseObject* parent=(BaseObject*)data->GetLink(CTT_OBJECT_LINK,doc,Obase);
    if (!parent) return NULL;
    
    Real maxSeg = data->GetReal(CTTSPOBJECT_MAXSEG,30.);
    LONG crntFrame = doc->GetTime().GetFrame(doc->GetFps());

    LONG delta = data->GetLong(CTTSPOBJECT_WINDOW,1);
    LONG splinePercentage = data->GetLong(CTT_SPLINE_PERCENTAGE,1);

    LONG strtFrame = crntFrame - delta;
    strtFrame = strtFrame<0?0:strtFrame;

    GeDynamicArray<BaseObject*> children;

    BaseObject* chld = NULL;
    LONG trck = 0;
    for (chld=parent->GetDownLast(); chld; chld=chld->GetPred()) {
        if (trck > strtFrame && trck<= crntFrame){
            children.Push((BaseObject*)chld->GetClone(COPYFLAGS_NO_HIERARCHY|COPYFLAGS_NO_ANIMATION|COPYFLAGS_NO_BITS,NULL));
        }
        trck++;
    }
    if (children.GetCount() < 2) {
        splineAtPoint.FreeArray();
        prvsFrame = crntFrame;
        return NULL;
    }
    
    if (children.GetCount() < delta && prvsFrame > crntFrame) {
        splineAtPoint.FreeArray();
    }

    SplineObject* ret = GetSpline(op, bt, doc, children, maxSeg, delta, splinePercentage);
    
    for (int k=0; k<children.GetCount(); k++){
        if (children[k]){
            BaseObject::Free(children[k]);
        }
    }
    prvsFrame = crntFrame;
    return ret;
}

SplineObject* SplinifyData::GetSpline(BaseObject* op, BaseThread* bt, BaseDocument *doc, GeDynamicArray<BaseObject*> &children, Real maxSeg, LONG delta, LONG splinePercentage){
    LONG child_cnt = children.GetCount();
    BaseContainer *data = op->GetDataInstance();
    
    LONG splineInterpolation = data->GetLong(SPLINEOBJECT_INTERPOLATION);

	StatusSetBar(0);
    StatusSetText("Collecting Points");
    GeDynamicArray<KDNode*> trees(child_cnt);
    GeDynamicArray<GeDynamicArray<Vector> > chldPoints(child_cnt);
    
    rng.Init(1244);
    LONG crntFrame = doc->GetTime().GetFrame(doc->GetFps());
    LONG shift = crntFrame - prvsFrame;
    shift = shift < 0? 0: shift;
    LONG startChild = prvsFrame == 0 || crntFrame == prvsFrame? 0 : child_cnt - shift - 1;
    startChild = startChild < 0? 0 : startChild;
    LONG maxPointCnt = 0;

    for (int k=splineAtPoint.GetCount() == 0 ? startChild: startChild+1; k < child_cnt; k++){
        Matrix ml;
        DoRecursion(op,children[k],chldPoints[k], ml);
        KDNode *kdTree;
        buildKDTree(chldPoints[k], &kdTree, rng);
        trees[k] = kdTree;
        if (maxPointCnt < chldPoints[k].GetCount()){
            maxPointCnt = chldPoints[k].GetCount();
        }
    }
    
    StatusSetBar(5);
    StatusSetText("Connecting Points");
    
    Real distMin = MAXREALr;
    Real distMax = 0.;
    
    for (LONG i = 0; i < splineAtPoint.GetCount(); i++){
        LONG shift = splineAtPoint[i].GetCount() - delta;
        if (shift > 0) {
            splineAtPoint[i].Shift(0, -shift);
            splineAtPoint[i].ReSize(splineAtPoint[i].GetCount() - shift);
        }
    }
    
    SplineObject* emptySpline = SplineObject::Alloc(0, SPLINETYPE_AKIMA);
    
    std::vector<SplinePair >splinePairs;

    Real avSplineSize = 0.0, avSplineLength = 0.0;

    if (splineAtPoint.GetCount() == 0){
        Random r;
        r.Init(43432);
        GeDynamicArray<LONG> validPoints(chldPoints[startChild].GetCount());
        validPoints.Fill(0,chldPoints[startChild].GetCount(),1);
        LONG fraction = chldPoints[startChild].GetCount()*splinePercentage/100;
        LONG cnt = 0, cycleCount = 0;
        while (cnt < fraction || cycleCount < 2*fraction) {
            cycleCount++;
            LONG indx = chldPoints[startChild].GetCount()*r.Get01();
            Vector queryPoint = chldPoints[startChild][indx];
            
            Real dist = -1.;
            LONG closestIndx = trees[startChild+1]->getNearestNeighbor(chldPoints[startChild+1], queryPoint, validPoints, dist, 0);
            if (closestIndx == -1){
                GePrint("error finding neighbor");
                continue;
            }
            
            if (dist > maxSeg || dist < 0.01) {
                continue;
            }

            if (dist > 0.0){
                validPoints[closestIndx] = 0;
                GeDynamicArray<Vector> rawSpline;
                rawSpline.Push(queryPoint);
                splineAtPoint.Push(rawSpline);
                cnt++;
            }
        }
    }

    GeDynamicArray<GeDynamicArray<LONG> > validPoints(splineAtPoint.GetCount());
    for (LONG k=0; k < splineAtPoint.GetCount(); k++){
        validPoints[k] = GeDynamicArray<LONG>(maxPointCnt);
        validPoints[k].Fill(0,maxPointCnt,1);//child_cnt - startChild
    }
    //////////
    AutoAlloc<SplineHelp> splineHelp;
    for (LONG i = 0; i < splineAtPoint.GetCount(); i++){
        GeDynamicArray<Vector> queryVector = splineAtPoint[i];
        for (int k=startChild; k < child_cnt-1; k++){

            Vector queryPoint = queryVector[queryVector.GetCount()-1];
            
            Real dist = -1.;
            LONG closestIndx = trees[k+1]->getNearestNeighbor(chldPoints[k+1], queryPoint, validPoints[k], dist, 0);
            if(closestIndx == -1){
                GePrint("error finding neighbor");
                continue;
            }
            
            distMin = distMin < dist ? distMin : dist;
            distMax = distMax > dist ? distMax : dist;
            
            if (dist > maxSeg || dist < 0.01) {
                continue;
            }
            validPoints[k][closestIndx] = 0;
            Vector clsst = chldPoints[k+1][closestIndx];
            if (splineAtPoint[i].Find(clsst) == NOTOK){
                splineAtPoint[i].Push(clsst);
            }
        }
        if (splineAtPoint[i].GetCount() == 0) continue;
        
        SplineObject	*spline=SplineObject::Alloc(splineAtPoint[i].GetCount(),SPLINETYPE_LINEAR);
        if (!spline) continue;
        
        spline->GetDataInstance()->SetBool(SPLINEOBJECT_CLOSED, FALSE);
        
        Vector *padr = spline->GetPointW();
        for (LONG l=0;l<splineAtPoint[i].GetCount();l++){
            padr[l] = splineAtPoint[i][l];
        }

        splineHelp->InitSpline(spline);
        Real splnLength = splineHelp->GetSplineLength();
        if (splnLength>0.0){
            spline->InsertUnder(emptySpline);
            splinePairs.push_back(SplinePair(spline, splineHelp->GetSplineLength()));
            avSplineLength += splnLength;
            avSplineSize += splineAtPoint[i].GetCount();
        }        

        if(i % 5 == 0){
            StatusSetBar(10 + (90*i)/splineAtPoint.GetCount());
            if (bt && bt->TestBreak()){
                break;
            }
        }
    }
    
//    std::sort(splinePairs.begin(), splinePairs.end(),comparator);
//    LONG limit = splinePairs.size()<80?splinePairs.size():80;
//    for (int s = 0; s < limit; s++){
//        std::cout<<splinePairs[s].second<<std::endl;
//        splinePairs[s].first->InsertUnder(emptySpline);
//    }

    String sizeAvg = splineAtPoint.GetCount() == 0? "Nan":RealToString(avSplineSize/splineAtPoint.GetCount());
    
    GePrint("pnts="+sizeAvg+" wndw="+LongToString(startChild)+"-" + LongToString(child_cnt)+" d="+RealToString(distMin)+":"+RealToString(distMax)+" sl="+RealToString(avSplineLength/avSplineSize));
    
    for (int k=0; k<child_cnt; k++){
        GeFree(trees[k]);
    }

	StatusClear();
    
    if (splinePairs.size() == 0) return NULL;
    
    ModelingCommandData mcd;
    mcd.doc = doc;
    mcd.op = emptySpline;
    
    if(!SendModelingCommand(MCOMMAND_JOIN, mcd)){
        return NULL;
    }
    
    return ToSpline(mcd.result->GetIndex(0L));
    
Error:
    for (int i = 0; i < children.GetCount(); i++){
        BaseObject::Free(children[i]);
    }
	return NULL;
}


Bool RegisterSplinify(void)
{
	return RegisterObjectPlugin(ID_SPLINIFYOBJECT,GeLoadString(IDS_SPLINIFY),OBJECT_GENERATOR|OBJECT_INPUT|OBJECT_ISSPLINE|OBJECT_CALL_ADDEXECUTION,SplinifyData::Alloc,"Octsplinify",AutoBitmap("tsp.tif"),0);
}

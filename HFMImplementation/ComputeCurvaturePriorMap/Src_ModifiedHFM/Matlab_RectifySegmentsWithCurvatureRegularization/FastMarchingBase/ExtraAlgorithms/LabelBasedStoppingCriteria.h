#ifndef LabelBasedStoppingCriteria_h
#define LabelBasedStoppingCriteria_h

#include "Base/HFMInterface.h"

template<typename T> struct LabelBasedStoppingCriteria :
HamiltonFastMarching<T>::ExtraAlgorithmInterface {
    typedef HamiltonFastMarching<T> HFM;
    typedef typename HFM::ExtraAlgorithmInterface Superclass;
    Redeclare7Types(HFM,IndexCRef,IndexType,ScalarType,Traits,Decision,HFMI,PointType)
    Redeclare1Constant(HFM,Dimension)
    
    virtual void Setup(HFMI*) override;
    virtual void Finally(HFMI*) override;
    virtual bool ImplementIn(HFM*) override;
    
    //------------------------------------------//
    std::map<IndexType,int,typename IndexType::LexicographicCompare> stopWhenLabelsAccepted;
    std::vector<IndexType> chosenTipsPerLabel;
    std::vector<double> labelsAtChosenTips;
    std::vector<double> valuesAtChosenTips;
    std::set<int> listLabels;
    //------------------------------------------//
    
    typedef LinearAlgebra::Point<ScalarType, Dimension+1> TaggedPointType;
//    typedef typename HFMI::SpecializationsDefault HFMIS;
//    typedef typename HFMIS::UnorientedIndexType UnorientedIndexType;
    
protected:
    virtual int PostProcess(IndexCRef) override;
    const HFM*  pFM = nullptr;
};



template<typename T> bool
LabelBasedStoppingCriteria<T>::ImplementIn(HFM * _pFM) {
    if(stopWhenLabelsAccepted.empty())
        return false;

    _pFM->extras.postProcess.push_back(this);
    pFM = _pFM;
    
    return true;
}

template<typename T> int
LabelBasedStoppingCriteria<T>::PostProcess(IndexCRef acceptedIndex){
    assert(pFM!=nullptr);
    
    const auto search=stopWhenLabelsAccepted.find(acceptedIndex);
    if(search!=stopWhenLabelsAccepted.end()){
        int label=search->second;
        if(listLabels.find(label)!= listLabels.end()){
            listLabels.erase(label);
            labelsAtChosenTips.push_back(1.0*label);
            chosenTipsPerLabel.push_back(acceptedIndex);
            valuesAtChosenTips.push_back(pFM->values(acceptedIndex));
        }
    }
    if(listLabels.empty())
        return Decision::kTerminate;
    
    return Decision::kAccept;
}

template<typename T> void
LabelBasedStoppingCriteria<T>::Setup(HFMI*that){
    auto & io = that->io;

    const auto & dom = that->pFM->dom;
    const auto & param = that->stencil.Param();
    
    // Stopping criteria
    if( io.HasField("stopWhenLabelsAccepted") ){
        const std::vector<TaggedPointType>& targets = io.template GetVector<TaggedPointType>("stopWhenLabelsAccepted");
        for(const auto & tp : targets) {
            PointType p;
            for(int k=0;k<Dimension;k++)
                p[k]=tp[k];
            int tag=(int)(tp[Dimension]);
            listLabels.insert(tag);
            IndexType index = dom.IndexFromPoint(param.ADim(p));
            if(dom.Periodize(index,index).IsValid())
                stopWhenLabelsAccepted.insert({index,tag});
            else
                Msg() << "Warning : target " << p << " is out of range.";
        }
    }
        
};

template<typename T> void
LabelBasedStoppingCriteria<T>::Finally(HFMI* that){
    
    if(stopWhenLabelsAccepted.empty())
        return;
    
    that->io.SetVector("labelsAtChosenTips",labelsAtChosenTips);
    that->io.SetVector("valuesAtChosenTips",valuesAtChosenTips);
    
    std::vector<PointType> chosenPoints;
    for(const IndexType& idx:chosenTipsPerLabel){
        PointType chosenPoint=this->pFM->dom.PointFromIndex(idx);
        chosenPoints.push_back(chosenPoint);
    }
    that->ExportGeodesics("",chosenPoints);
    
    std::vector<PointType> chosenOutputPoints;
    for(const PointType& pt:chosenPoints){
        PointType newPt;
        for(int k=0;k<Dimension;k++)
            newPt[k]=pt[k]-0.5;
        chosenOutputPoints.push_back(that->stencil.Param().ReDim(newPt));
    }
    
    that->io.SetVector("chosenTipsPerLabel", chosenOutputPoints);
    that->io.SetVector("chosenTipsPerLabelIndex",chosenPoints);
};
#endif

#ifndef LabelRetrieval_h
#define LabelRetrieval_h


template<typename T> struct LabelRetrieval :
HamiltonFastMarching<T>::ExtraAlgorithmInterface {
    typedef HamiltonFastMarching<T> HFM;
    typedef typename HFM::ExtraAlgorithmInterface Superclass;
    Redeclare13Types(HFM,IndexCRef,IndexType,ScalarType,Traits,HFMI,PointType,
					 DiscreteType,OffsetCRef,VectorType,DiscreteFlowType,RecomputeType,
					 DiscreteFlowElement,Decision)
    Redeclare1Constant(HFM,Dimension)
     
    typedef LinearAlgebra::Point<ScalarType, Dimension+1> TaggedPointType;
    typename HFM::template Array<ScalarType,Dimension> labelMap;
    
    virtual void Setup(HFMI*) override;
    virtual void Finally(HFMI*) override;
    virtual bool ImplementIn(HFM*) override;
protected:
    std::vector<TaggedPointType> taggedPoints;
    virtual int PostProcess(IndexCRef) override;
    const HFM * pFM=nullptr;
};

template<typename T> void LabelRetrieval<T>::
Setup(HFMI * that){
    auto & io = that->io;
    if(io.HasField("labelMap"))
        labelMap=io.template GetArray<ScalarType, Dimension>("labelMap");
    else
        labelMap.clear();
}

template<typename T> void LabelRetrieval<T>::
Finally(HFMI* that){
    auto & io = that->io;
    io.SetVector("taggedPoints",taggedPoints);
}

template<typename T> bool LabelRetrieval<T>::
ImplementIn(HFM*_pFM){
    if(labelMap.empty())
        return false;
    
    pFM=_pFM;
    _pFM->extras.postProcess.push_back(this);
    return true;
}

template<typename T> int LabelRetrieval<T>::
PostProcess(IndexCRef index){
    assert(pFM!=nullptr);
    assert(!labelMap.empty());
    ScalarType label = labelMap(index);
    if(label>0.5){
        TaggedPointType taggedPt;
        for(int k=0;k<Dimension;k++)
            taggedPt[k]=index[k];
        taggedPt[Dimension]=label;
        taggedPoints.push_back(taggedPt);
    }
    return 0;
}

#endif /* LabelRetrieval_h */

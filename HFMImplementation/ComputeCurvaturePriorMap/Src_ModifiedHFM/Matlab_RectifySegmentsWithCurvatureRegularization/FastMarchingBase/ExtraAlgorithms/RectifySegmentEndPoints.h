#ifndef RectifySegmentEndPoints_h
#define RectifySegmentEndPoints_h

template<typename T> struct RectifySegmentEndPoints :
HamiltonFastMarching<T>::ExtraAlgorithmInterface {
    typedef HamiltonFastMarching<T> HFM;
    typedef typename HFM::ExtraAlgorithmInterface Superclass;
    Redeclare13Types(HFM,IndexCRef,IndexType,ScalarType,Traits,HFMI,PointType,
					 DiscreteType,OffsetCRef,VectorType,DiscreteFlowType,RecomputeType,
					 DiscreteFlowElement,Decision)
    Redeclare1Constant(HFM,Dimension)
    
    ScalarType admEuclidLength=0.0;
    IndexType rectifiedIndex = IndexType::Constant(-1);
    std::unique_ptr<typename HFM::template DataSource<VectorType> > euclideanScale;
    typename HFM::template Array<ScalarType,Dimension> euclideanLengths;
    
    virtual void Setup(HFMI*) override;
    virtual void Finally(HFMI*) override;
    virtual bool ImplementIn(HFM*) override;
protected:
    bool exportRectifiedPoint=false;
    virtual int PostProcessWithRecompute(IndexCRef, const RecomputeType &, const DiscreteFlowType &) override;
    const HFM * pFM=nullptr;
};

template<typename T> void RectifySegmentEndPoints<T>::
Setup(HFMI * that){
    auto & io = that->io;
    const auto & dom = that->pFM->dom;
    if(io.HasField("exportRectifiedPoint"))
        exportRectifiedPoint=io.template Get<ScalarType>("exportRectifiedPoint")>0.5;
    else
        exportRectifiedPoint=false;
    
    if(io.HasField("euclideanScale"))
        euclideanScale = that->template GetField<VectorType>("euclideanScale");
    
    assert(io.HasField("admEuclidLength"));
    if(io.HasField("admEuclidLength"))
        admEuclidLength = io.template Get<ScalarType>("admEuclidLength");
    
}

template<typename T> void RectifySegmentEndPoints<T>::
Finally(HFMI* that){
    auto & io = that->io;
    if(exportRectifiedPoint){
        io.template Set<PointType>("rectifiedIndex",PointType::CastCoordinates(rectifiedIndex));
        that->ExportGeodesics("rectified",{that->pFM->dom.PointFromIndex(rectifiedIndex)});
    }
}

template<typename T> bool RectifySegmentEndPoints<T>::
ImplementIn(HFM*_pFM){
    if(!exportRectifiedPoint)
        return false;
    
    pFM=_pFM;
    _pFM->extras.postProcessWithRecompute.push_back(this);
    euclideanLengths.dims=pFM->values.dims;
    euclideanLengths.resize(euclideanLengths.dims.Product(),-1);
    return true;
}

template<typename T> int RectifySegmentEndPoints<T>::
PostProcessWithRecompute(IndexCRef index, const RecomputeType &, const DiscreteFlowType & flow){
    assert(pFM!=nullptr);
    
    VectorType flowSum=VectorType::Constant(0);
    ScalarType lengthSum=0, weightSum=0;
    for(const DiscreteFlowElement & fl : flow){
        if(fl.weight==0) continue;
        weightSum+=fl.weight;
        flowSum+=fl.weight * VectorType::CastCoordinates(fl.offset);
        IndexType neigh = index;
        for(int i=0; i<Dimension; ++i)
            neigh[i]+=fl.offset[i];
        const auto transform = pFM->dom.Periodize(neigh,index);
        assert(transform.IsValid());
        (void)transform;
        lengthSum+=fl.weight * euclideanLengths(neigh);
    }
    if(weightSum>0){
        flowSum/=weightSum;
        lengthSum/=weightSum;
        if(euclideanScale==nullptr){
            const VectorType scale=VectorType::Constant(1.0);
            for(int i=0; i<Dimension; ++i)
                flowSum[i]*=scale[i];
            const ScalarType length = flowSum.Norm() + lengthSum;
            euclideanLengths(index) = length;
            if(length>=admEuclidLength){
                rectifiedIndex=index;
                return Decision::kTerminate;
            }
        }
        else{
            const VectorType scale = (*euclideanScale)(index);
            for(int i=0; i<Dimension; ++i)
                flowSum[i]*=scale[i];
            const ScalarType length = flowSum.Norm() + lengthSum;
            euclideanLengths(index) = length;
            
            if(length>=admEuclidLength){
                rectifiedIndex=index;
                return Decision::kTerminate;
            }
        }
    }
    else{
        euclideanLengths(index) = 0;
    }
    
    return 0;
}

#endif /* RectifySegmentEndPoints_h */

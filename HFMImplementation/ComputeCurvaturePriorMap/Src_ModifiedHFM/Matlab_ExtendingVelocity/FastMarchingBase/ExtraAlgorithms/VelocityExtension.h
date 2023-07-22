#ifndef VelocityExtension_h
#define VelocityExtension_h

template<typename T> struct VelocityExtension :
HamiltonFastMarching<T>::ExtraAlgorithmInterface {
    typedef HamiltonFastMarching<T> HFM;
    typedef typename HFM::ExtraAlgorithmInterface Superclass;
    Redeclare15Types(HFM,IndexCRef,IndexType,ScalarType,Traits,HFMI,PointType,
					 DiscreteType,ShortType,OffsetCRef,VectorType,DiscreteFlowType,
					 RecomputeType,DiscreteFlowElement,Decision,IndexDiff)
    Redeclare1Constant(HFM,Dimension)
    
    typename HFM::template Array<ScalarType,Dimension> extendedVelocity;
    
    virtual void Setup(HFMI*) override;
    virtual void Finally(HFMI*) override;
    virtual bool ImplementIn(HFM*) override;
protected:
    virtual int PostProcessWithRecompute(IndexCRef, const RecomputeType &, const DiscreteFlowType &) override;
    const HFM * pFM=nullptr;
};

template<typename T> void VelocityExtension<T>::
Setup(HFMI * that){
    auto & io = that->io;
    if(!(io.HasField("initialVelocity")))
        return;
    extendedVelocity.dims = that->stencil.dims;
    extendedVelocity.resize(that->stencil.dims.Product(),0);
    
    if(io.HasField("initialVelocity")){
        const auto initialVelocity = io.template GetVector<ScalarType>("initialVelocity");
        const auto seeds = io.template GetVector<PointType>("seeds");
        if(initialVelocity.size()!=seeds.size()){
            ExceptionMacro("Error : Number of seeds (" <<seeds.size()<<
                           ") is distinct from the number of initialVelocity(" <<initialVelocity.size() << ").\n");}
        for(size_t i=0; i<seeds.size(); ++i){
            IndexType index = that->pFM->dom.IndexFromPoint(that->stencil.Param().ADim(seeds[i]));
            if(!that->pFM->dom.Periodize(index,index).IsValid())
                ExceptionMacro("Error : seed " << seeds[i] << " is out of range.\n");
            extendedVelocity(index) = (ScalarType)initialVelocity[i];
        }
    }
    
}

template<typename T> bool VelocityExtension<T>::
ImplementIn(HFM*_pFM){
    if(extendedVelocity.empty())
        return false;
    
    _pFM->extras.postProcessWithRecompute.push_back(this);
    pFM=_pFM;
    assert(extendedVelocity.dims==pFM->values.dims);
    return true;
}

template<typename T> void VelocityExtension<T>::
Finally(HFMI*that){
    auto & io = that->io;
    if(io.template Get<ScalarType>("exportExtendedVelocity",1))
        io.SetArray("extendedVelocity",extendedVelocity.template Cast<ScalarType>());
}

template<typename T> int VelocityExtension<T>::
PostProcessWithRecompute(IndexCRef index, const RecomputeType &, const DiscreteFlowType & flow) {
    
    ScalarType & currentVelocity = extendedVelocity(index);
    ScalarType maxWeight=0;
    for(const DiscreteFlowElement & fl : flow){
        if(fl.weight<=maxWeight)
            continue;
        maxWeight=fl.weight;
        IndexType neigh = index + IndexDiff::CastCoordinates(fl.offset);
        const auto transform = pFM->dom.Periodize(neigh,index);
        assert(transform.IsValid());
        (void)transform;
        currentVelocity = extendedVelocity(neigh);
    }
    
    return 0;
}

#endif /* VoronoiDiagram_h */

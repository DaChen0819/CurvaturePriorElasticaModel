#ifndef RegionalCurvature2_h
#define RegionalCurvature2_h


#include "Specializations/Curvature2.h"
#include "Experimental/RegionalCurvature2.h"

//// --------------- R2S1NonShared Traits ---------------
//
//struct TraitsR2S1NonShared : TraitsBase<3> {
//    typedef EulerianDifference<OffsetType, ScalarType, 0> DifferenceType;
//    constexpr static std::array<Boundary, Dimension> boundaryConditions =
//    {{Boundary::Closed, Boundary::Closed, Boundary::Periodic}};
//	using DomainType = PeriodicGrid<TraitsR2S1NonShared>;
//    // Stencils actually depend on all coordinates (no multiplier). This is to get proper domain parametrization.
//    static const DiscreteType nStencilDependencies=1;
//    constexpr static std::array<DiscreteType, nStencilDependencies> stencilDependencies = {{2}};
//};
//
//constexpr decltype(TraitsR2S1NonShared::boundaryConditions) TraitsR2S1NonShared::boundaryConditions;
//constexpr decltype(TraitsR2S1NonShared::stencilDependencies) TraitsR2S1NonShared::stencilDependencies;

// ----------- 2D Reeds-Shepp Regional model ---------
struct TraitsReedsSheppRegion2 : TraitsR2S1NonShared {
    typedef EulerianStencil<DifferenceType,6+1> StencilType;
};

struct StencilReedsSheppRegion2 final
: HamiltonFastMarching<TraitsReedsSheppRegion2>::StencilDataType {
    typedef HamiltonFastMarching<TraitsReedsSheppRegion2> HFM;
    typedef HFM::StencilDataType Superclass;
    HFM::ParamDefault param;
    ScalarType eps=0.1;
    typedef Traits::DataSource<ScalarType> ScalarFieldType;
    std::unique_ptr<ScalarFieldType> pSpeed, pXi, pTheta;
	typedef Traits::BasisReduction<3> ReductionType;
	Voronoi1Vec<ReductionType> reduc;
	
    virtual void SetStencil(IndexCRef index, StencilType & stencil) override {
		// This implementation could be improved, probably, by using a single tensor
		// (One less offset, and probably smaller offsets.)
        assert(pSpeed); assert(pXi);
        const ScalarType
        speed=(*pSpeed)(index),
		xi=(*pXi)(index),
        gS=param.gridScale,
        tS=param.dependScale;
		const ScalarType theta = pTheta ? (*pTheta)(index) : index[2]*tS;
        const ScalarType c = cos(theta), s=sin(theta);
        const VectorType v{c/gS,s/gS,0.0};
        
        auto & symmetric = stencil.symmetric[0];
        reduc(&symmetric[0], speed*v);
        
        symmetric[6].offset = OffsetType{0,0,1};
        symmetric[6].baseWeight = square(speed/(xi*tS));
    }
    virtual const ParamInterface & Param() const override {return param;}
    virtual void Setup(HFMI *that) override {
        Superclass::Setup(that);
		auto & io = that->io;
        reduc.eps=io.template Get<ScalarType>("eps",reduc.eps);
        typedef typename HFMI::template DataSource_Inverse<ScalarType> SourceInvType;
        if(io.HasField("speed")) pSpeed = that->template GetField<ScalarType>("speed",false);
        else pSpeed = std::unique_ptr<SourceInvType>(new SourceInvType(that->template GetField<ScalarType>("cost",false) ) );
        pXi = that->GetField<ScalarType>("xi",false);
		if(io.HasField("theta")){
            pTheta = that->GetField<ScalarType>("theta",false);
        }
        param.Setup(io,2*mathPi/dims.back());}
};

// ----------- 2D Reeds-Shepp Forward Extended model ---------

struct TraitsReedsSheppForwardRegion2 : TraitsR2S1NonShared {
    typedef EulerianStencil<DifferenceType,1,6> StencilType;
};

struct StencilReedsSheppForwardRegion2 final
: HamiltonFastMarching<TraitsReedsSheppForwardRegion2>::StencilDataType {
    typedef HamiltonFastMarching<TraitsReedsSheppForwardRegion2> HFM;
    typedef HFM::StencilDataType Superclass;
    HFM::ParamDefault param;
    ScalarType eps=0.1;
    typedef Traits::DataSource<ScalarType> ScalarFieldType;
    std::unique_ptr<ScalarFieldType> pSpeed, pXi, pTheta;
	typedef Traits::BasisReduction<3> ReductionType;
	Voronoi1Vec<ReductionType> reduc;
    virtual void SetStencil(IndexCRef index, StencilType & stencil) override {
        assert(pSpeed); assert(pXi); assert(pKappa);
        const ScalarType
        speed=(*pSpeed)(index),
		xi=(*pXi)(index),
        gS=param.gridScale,
        tS=param.dependScale;
		const ScalarType theta = pTheta ? (*pTheta)(index) : index[2]*tS;
        const ScalarType c = cos(theta), s=sin(theta);
        const VectorType v{c/gS,s/gS,0.0};
        
        auto & forward = stencil.forward[0];
        reduc(&forward[0], speed*v);
        
        auto & symmetric = stencil.symmetric[0];
        symmetric[0].offset = OffsetType{0,0,1};
        symmetric[0].baseWeight = square(speed/(xi*tS));
    }
    virtual const ParamInterface & Param() const override {return param;}
    virtual void Setup(HFMI *that) override {
        Superclass::Setup(that);
		auto & io = that->io;
        reduc.eps=io.template Get<ScalarType>("eps",reduc.eps);
        typedef typename HFMI::template DataSource_Inverse<ScalarType> SourceInvType;
        if(io.HasField("speed")) pSpeed = that->template GetField<ScalarType>("speed",false);
        else pSpeed = std::unique_ptr<SourceInvType>(new SourceInvType(that->template GetField<ScalarType>("cost",false) ) );
        pXi = that->GetField<ScalarType>("xi");
		if(io.HasField("theta")) {
            pTheta = that->GetField<ScalarType>("theta",false);
        }
        param.Setup(io,2*mathPi/dims.back());}
};


// ---------- 2D Dubins Extended model -------------

struct TraitsDubinsRegion2 : TraitsR2S1NonShared {
    typedef EulerianStencil<DifferenceType,0,6,2> StencilType;
};

struct StencilDubinsRegion2 final
: HamiltonFastMarching<TraitsDubinsRegion2>::StencilDataType {
    typedef HamiltonFastMarching<TraitsDubinsRegion2> HFM;
    typedef HamiltonFastMarching<TraitsDubinsRegion2>::StencilDataType Superclass;
    HFM::ParamDefault param;
    ScalarType eps=0.1;
    typedef Traits::DataSource<ScalarType> ScalarFieldType;
    std::unique_ptr<ScalarFieldType> pSpeed, pXi, pTheta;
	typedef Traits::BasisReduction<3> ReductionType;
	Voronoi1Vec<ReductionType> reduc;
	
    virtual void SetStencil(IndexCRef index, StencilType & stencil) override {
        assert(pSpeed); assert(pXi); assert(pKappa);
        const ScalarType
        speed=(*pSpeed)(index),
		xi=(*pXi)(index),
        gS=param.gridScale,
        tS=param.dependScale;
		const ScalarType theta = pTheta ? (*pTheta)(index) : index[2]*tS;
        const ScalarType c = cos(theta), s=sin(theta);
        const VectorType
        vL{c/gS,s/gS,(1./xi)/tS},
        vR{c/gS,s/gS,(-1./xi)/tS};
        
        reduc(&stencil.forward[0][0], speed*vL);
        reduc(&stencil.forward[1][0], speed*vR);
	}
    virtual const ParamInterface & Param() const override {return param;}
    virtual void Setup(HFMI *that) override {
        Superclass::Setup(that);
		auto & io = that->io;
		reduc.eps=io.template Get<ScalarType>("eps",reduc.eps);
        typedef typename HFMI::template DataSource_Inverse<ScalarType> SourceInvType;
        if(io.HasField("speed")) pSpeed = that->template GetField<ScalarType>("speed",false);
        else pSpeed = std::unique_ptr<SourceInvType>(new SourceInvType(that->template GetField<ScalarType>("cost",false) ) );
        pXi = that->GetField<ScalarType>("xi");
		if(io.HasField("theta")){
            pTheta = that->GetField<ScalarType>("theta",false);
        }
        param.Setup(io,2*mathPi/dims.back());}
};


// ----------- 2D Euler-Mumford Elastica Extended model -----------

template<int nFejer>
struct TraitsElasticaRegion2 : TraitsR2S1NonShared {
	using DifferenceType = typename TraitsR2S1NonShared::DifferenceType;
    using StencilType = EulerianStencil<DifferenceType,0,6*nFejer>;
};

template<int nFejer>
struct StencilElasticaRegion2 final
: HamiltonFastMarching<TraitsElasticaRegion2<nFejer> >::StencilDataType {
    typedef HamiltonFastMarching<TraitsElasticaRegion2<nFejer> > HFM;
    typedef typename HFM::StencilDataType Superclass;
    Redeclare7Types(Superclass,Traits,ScalarType,IndexCRef,
					VectorType,StencilType,ParamInterface,HFMI)
    Redeclare1Constant(Traits,mathPi)
    typename HFM::ParamDefault param;
    ScalarType eps=0.1;
    typedef typename Traits::template DataSource<ScalarType> ScalarFieldType;
    std::unique_ptr<ScalarFieldType> pSpeed, pXi, pTheta;
	typedef typename Traits::template BasisReduction<3> ReductionType;
	Voronoi1Vec<ReductionType> reduc;
	
    virtual void SetStencil(IndexCRef index, StencilType & stencil) override {
        assert(pSpeed); assert(pXi); assert(pKappa);
        const ScalarType
        speed=(*pSpeed)(index),
		xi=(*pXi)(index),
        gS=param.gridScale,
        tS=param.dependScale;
		const ScalarType theta = pTheta ? (*pTheta)(index) : index[2]*tS;
        const ScalarType cT = cos(theta), sT=sin(theta);
        
        auto & forward = stencil.forward[0];
        for(int l=0; l<nFejer; ++l){
            const ScalarType phi = mathPi*(l+0.5)/nFejer;
            const ScalarType cP = cos(phi), sP=sin(phi);
            const VectorType v{sP*cT/gS,sP*sT/gS,(cP/xi)/tS};

            reduc(&forward[6*l], v*speed);
            for(int i=0; i<6; ++i) forward[6*l+i].baseWeight*=StencilElastica2<nFejer>::fejerWeights[l];
        }
    }
    virtual const ParamInterface & Param() const override {return param;}
    virtual void Setup(HFMI *that) override {
        Superclass::Setup(that);
		auto & io = that->io;
        reduc.eps=io.template Get<ScalarType>("eps",reduc.eps);
        typedef typename HFMI::template DataSource_Inverse<ScalarType> SourceInvType;
        if(io.HasField("speed")) pSpeed = that->template GetField<ScalarType>("speed",false);
        else pSpeed = std::unique_ptr<SourceInvType>(new SourceInvType(that->template GetField<ScalarType>("cost",false) ) );
        pXi = that->template GetField<ScalarType>("xi");
		if(io.HasField("theta")){
            pTheta = that->template GetField<ScalarType>("theta",false);
        }
        param.Setup(io,2*mathPi/this->dims.back());}
};

#endif /* RegionalCurvature2_h */

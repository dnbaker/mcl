#include "blaze/Math.h"

namespace bmcl {
using blaze::unchecked;

struct MCLSettings {
    unsigned matrix_expand_ = 2;
    double inflate_ = 2.;
    double mult_    = 1.;
    long unsigned niter_ = 1000;
    double selfloop_inc_ = 1;
    //double threshold_ = 1e-9; // Values < this are set to 0 at each iteration
};

template<typename MT>
void mcl(blaze::Matrix<MT, blaze::rowMajor> &matrix, MCLSettings settings) {
    // Store the transpose of the matrix
    // Such that each column corresponds to
    // diffusion from one's self to others
    auto &mref = (*matrix);
    using FT = blaze::ElementType_t<MT>;
    // Add diagonal, if selected
    if(settings.selfloop_inc_ != 0.) {
        if constexpr(blaze::IsDenseMatrix_v<MT>) {
            band(mref, 0) += settings.selfloop_inc_;
        } else {
            for(unsigned i = 0; i < (mref).rows(); ++i) {
                auto r = row(mref, i, unchecked);
                auto it = r.find(i);
                if(it == r.end())
                    r.set(i, settings.selfloop_inc_);
                else
                    it->value() += settings.selfloop_inc_;
            }
        }
    }
    blaze::DynamicVector<FT, blaze::rowVector> isums;
    auto normalize = [&mref,&isums]() {
        isums = 1. / blaze::sum<blaze::columnwise>(mref);
        mref %= blaze::expand(isums, (mref).rows());
    };
    normalize();
    size_t niter = 0;
    for(;;) {
        // inflate
        (mref) = pow(mref, settings.inflate_);
        normalize();
        // expand
        mref *= mref;
        
        if(++niter == settings.niter_) break;
    }
    // TODO: emit clusters and terminate
}

} // namespace bmcl

#include "blaze/Math.h"

namespace bmcl {
using blaze::unchecked;

struct MCLSettings {
    unsigned matrix_expand_ = 2;
    double inflate_ = 2.;
    double mult_    = 1.;
    long unsigned niter_ = 1000;
    double selfloop_inc_ = 1;
    double threshold_ = 1e-9; // Values < this are set to 0 at each iteration
};

template<typename MT, bool SO>
void mcl(blaze::Matrix<MT, SO> &matrix, MCLSettings settings) {
    using FT = blaze::ElementType_t<MT>;
    // Add diagonal, if selected
    if(settings.selfloop_inc_ != 0.) {
        if constexpr(blaze::IsDenseMatrix_v<MT>) {
            band(~matrix, 0 += settings.selfloop_inc_);
        } else {
            OMP_PFOR
            for(unsigned i = 0; i < (~matrix).rows(); ++i) {
                auto r = row(~matrix, i, unchecked);
                auto it = r.find(i);
                if(it == r.end())
                    r.set(i, settings.selfloop_inc_);
                else
                    it->value() += settings.selfloop_inc_;
            }
        }
    }
    blaze::DenseMatrix<FT> isums;
    auto normalize = [&matrix,&isums]() {
        isums = 1. / blaze::sum<blaze::rowwise>(~matrix);
        ~matrix %= blaze::expand(isums, (~matrix).rows());
    };
    normalize();
    size_t niter = 0;
    for(;;) {
        // inflate
        (~matrix) = pow(~matrix, settings.inflate_);
        normalize();
        // expand
        ~matrix *= ~matrix;
        if(settings.threshold_ != 0.) {
            // prune small values
           ~matrix = clamp(~matrix, settings.threshold_, std::numeric_limits<FT>::max());
        }
        if(++iter == settings.niter_) break;
        if(max(pow(~matrix, 2.) - ~matrix) - min(pow(~matrix, 2.) - ~matrix) <= settings.tol_)
            break;
    }
    // TODO: emit clusters and terminate
}

} // namespace bmcl

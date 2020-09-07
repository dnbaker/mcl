#include "blaze/Math.h"
#include "./macros.h"
#include <algorithm>

namespace bmcl {
using blaze::unchecked;

struct MCLSettings {
    unsigned matrix_expand_ = 2;
    double inflate_ = 2.;
    double mult_    = 1.;
    long unsigned niter_ = 1000;
    double selfloop_inc_ = 1;
    double threshold_ = 1e-7; // Values < this are set to 0 at each iteration
};

template<typename MT>
void prune(blaze::DenseMatrix<MT, blaze::rowMajor> &matrix, MCLSettings settings) {
    if(settings.threshold_ <= 0.) return;
    using FT = blaze::ElementType_t<MT>;
    auto fn = [t=settings.threshold_](auto x) {if(x <= t) x = FT(0.); return x;};
    *matrix = map(*matrix, fn);
}

template<typename MT>
void prune(blaze::SparseMatrix<MT, blaze::rowMajor> &matrix, MCLSettings settings) {
    if(settings.threshold_ <= 0.) return;
    auto fn = [t=settings.threshold_](auto x) {return x < t;};
    (*matrix).erase(fn);
}

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
        prune(mref, settings);
        normalize();
        // expand
        switch(settings.matrix_expand_) {
            case 4: mref *= mref; [[fallthrough]];
            case 2: mref *= mref; break;
            case 3: mref = mref * mref * mref; break;
            default: 
                std::fprintf(stderr, "Expansion values permitted: 2, 3, 4.\n");
                throw std::invalid_argument("matrix expand should be 2, 3, or 4");
        }
        
        if(++niter == settings.niter_) break;
    }
}

template<typename MT, typename IT=uint32_t>
auto get_clusters(const blaze::Matrix<MT, blaze::rowMajor> &matrix, MCLSettings settings) {
    // This effectively extracts out the non-empty rows
    // and returns them as sparse vectors
    // The rows correspond to cluster identities (contained in centers)
    // with the lists of items assigned in the sparse vectors of assignments
    auto &m = *matrix;
    blaze::DynamicVector<IT> centers(m.rows());
    size_t ncenters = 0;
    OMP_PFOR
    for(size_t i = 0; i < m.rows(); ++i) {
        auto r = row(m, i, unchecked);
        if(max(r) > 0.) {
            OMP_CRITICAL
            {
                centers[ncenters++] = i;
            }
        }
    }
    centers.resize(ncenters);
    std::sort(centers.data(), centers.data() + ncenters);
    std::vector<blaze::CompressedVector<blaze::ElementType_t<MT>>> assignments(centers.size());
    OMP_PFOR
    for(auto i = 0u; i < ncenters; ++i) {
        auto &asn = assignments[i];
        asn = row(m, i, unchecked);
    }
    return std::make_pair(centers, assignments);
}

} // namespace bmcl

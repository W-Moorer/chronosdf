#pragma once

#include <cstddef>
#include <vector>

#include "SdfNarrowphase.h"

namespace chrono_sdf_contact {

class SdfClustering {
  public:
    explicit SdfClustering(std::size_t max_contacts = 64, double cell_size = 0.05);

    void SetMaxContacts(std::size_t max_contacts);
    void SetCellSize(double cell_size);
    void Apply(std::vector<GeneratedContact>& contacts) const;

  private:
    std::size_t max_contacts_;
    double cell_size_;
};

}  // namespace chrono_sdf_contact

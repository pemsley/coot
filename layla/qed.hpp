//
//  Copyright (c) 2009-2017, Novartis Institutes for BioMedical Research Inc.,
//                2024, Jakub Smulski, Global Phasing Ltd.
//  All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//     * Neither the name of Novartis Institutes for BioMedical Research Inc.
//       nor the names of its contributors may be used to endorse or promote
//       products derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#ifndef LAYLA_QED_HPP
#define LAYLA_QED_HPP
#include <memory>
#include <vector>
#include <optional>
#include <rdkit/GraphMol/RWMol.h>

namespace coot::layla::RDKit {

class QED {
    public:
    struct QEDproperties {
        double MW;
        double ALOGP;
        double HBA;
        double HBD;
        double PSA;
        double ROTB;
        double AROM;
        double ALERTS;
    };

    enum class QEDPropName : std::size_t {
        MW = 0,
        ALOGP,
        HBA,
        HBD,
        PSA,
        ROTB,
        AROM,
        ALERTS
    };

    struct ADSparameter {
        double A;
        double B;
        double C;
        double D;
        double E;
        double F;
        double DMAX;
    };

    private:
    static const QEDproperties WEIGHT_MAX;
    static const QEDproperties WEIGHT_MEAN;
    static const QEDproperties WEIGHT_NONE;

    static const std::unique_ptr<const ::RDKit::ROMol> AliphaticRings;

    static const std::vector<std::unique_ptr<const ::RDKit::ROMol>> Acceptors;
    static const std::vector<std::unique_ptr<const ::RDKit::ROMol>> StructuralAlerts;
    /// Indexed via QEDPropName enum
    static const std::vector<ADSparameter> adsParameters;

    public:

   class QED_and_ads_t {
   public:
      double qed_score;
      double ads_mw;
      double ads_alogp;
      double ads_hba;
      double ads_hbd;
      double ads_psa;
      double ads_rotb;
      double ads_arom;
      double ads_alerts;

   };
   
    /// Calculates the properties that are required to calculate the QED descriptor.
    static QEDproperties properties(const ::RDKit::ROMol& mol);

   /// Calculate the weighted sum of ADS mapped properties
    // @setDescriptorVersion(version='1.1.0')
    static QED_and_ads_t qed(const ::RDKit::ROMol& mol,
                             std::optional<QEDproperties> qedProperties = std::nullopt,
                             QEDproperties w = WEIGHT_MEAN);

    /// ADS function
    static double ads(double x, const ADSparameter& p) noexcept;

    /// Calculates the QED descriptor using maximal descriptor weights
    inline static double weights_max(const ::RDKit::ROMol& mol) {
       return qed(mol, WEIGHT_MAX).qed_score;
    }

    /// Calculates the QED descriptor using average descriptor weights.
    inline static double weights_mean(const ::RDKit::ROMol& mol) {
       return qed(mol, WEIGHT_MEAN).qed_score;
    }

    /// Calculates the QED descriptor using unit weights.
    inline static double weights_none(const ::RDKit::ROMol& mol) {
       return qed(mol, WEIGHT_NONE).qed_score;
    }

    /// Calculates the QED descriptor using average descriptor weights.
    inline static double default_(const ::RDKit::ROMol& mol) {
        return weights_mean(mol);
    }
};


}

#endif

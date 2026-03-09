# Derive SCDM parameters (mu, sigma) from projectability fitting (Vitale et al. arXiv:1909.00433)

import xml.etree.ElementTree as ET
import numpy as np
from scipy.optimize import curve_fit
from scipy.special import erfc


def scdm_fit(
        structure_name,
        atomic_proj_xml=None
):
    """
    Parse projwfc.x output and return fitted SCDM parameters.
    structure_name : str
        Base name of the structure (used to locate default xml path).
    atomic_proj_xml : str, optional
        Explicit path to atomic_proj.xml. If None, defaults to
        '{structure_name}/{structure_name}.save/atomic_proj.xml'.
    Returns dict with keys: scdm_mu, scdm_sigma, mu_fit, sigma_fit
    """
    def parse_atomic_proj_xml(xml_path):
        """
        Returns arrays of eigenvalues (eV) and projectabilities for all (n,k) states.
        """
        tree = ET.parse(xml_path)
        root = tree.getroot()
        header = root.find('.//HEADER')
        nbnd = int(header.get('NUMBER_OF_BANDS'))
        eigenstates = root.find('.//EIGENSTATES')
        if eigenstates is None:
            raise ValueError(f"Could not find EIGENSTATES in {xml_path}")
        e_elems = list(eigenstates.iter('E'))
        projs_elems = list(eigenstates.iter('PROJS'))
        if len(e_elems) != len(projs_elems):
            raise ValueError(f"Mismatch: {len(e_elems)} E blocks vs {len(projs_elems)} PROJS blocks")
        eigenvalues = []
        projectabilities = []
        for e_elem, projs_elem in zip(e_elems, projs_elems):
            ek = np.array(e_elem.text.split(), dtype=float)
            eigenvalues.extend(ek)
            # p_nk = sum_i |<psi_nk|phi_i>|^2
            proj_k = np.zeros(nbnd)
            for wfc_elem in projs_elem.iter('ATOMIC_WFC'):
                vals = np.array(wfc_elem.text.split(), dtype=float)
                re_parts = vals[0::2]
                im_parts = vals[1::2]
                proj_k += re_parts**2 + im_parts**2
            projectabilities.extend(proj_k)
        eigenvalues = np.array(eigenvalues)
        projectabilities = np.array(projectabilities)
        ry_to_ev = 13.605693123
        eigenvalues *= ry_to_ev
        return eigenvalues, projectabilities

    def erfc_model(e, mu, sigma):
        return 0.5 * erfc((e - mu) / sigma)

    def fit_scdm_parameters(eigenvalues, projectabilities):
        """
        Fit erfc to (eigenvalue, projectability) data.
        Returns (scdm_mu, scdm_sigma, mu_fit, sigma_fit).
        """
        mu0 = np.median(eigenvalues)
        sigma0 = np.std(eigenvalues) / 4
        # f(e; mu, sigma) = 0.5 * erfc(-(mu - e) / sigma)
        popt, _ = curve_fit(
            erfc_model, eigenvalues, projectabilities,
            p0=[mu0, sigma0],
            bounds=([eigenvalues.min(), 0.01], [eigenvalues.max(), eigenvalues.max() - eigenvalues.min()]),
            maxfev=10000
        )
        mu_fit, sigma_fit = popt
        # scdm_mu = mu_fit - 3*sigma_fit, scdm_sigma = sigma_fit
        scdm_mu = mu_fit - 3.0 * sigma_fit
        scdm_sigma = sigma_fit
        return scdm_mu, scdm_sigma, mu_fit, sigma_fit

    # Args
    if atomic_proj_xml is None:
        atomic_proj_xml = f"{structure_name}/{structure_name}.save/atomic_proj.xml"
    # Parse and fit
    print(f"Parsing projectabilities from {atomic_proj_xml}...")
    eigenvalues, projectabilities = parse_atomic_proj_xml(atomic_proj_xml)
    print(f"  Collected {len(eigenvalues)} (eigenvalue, projectability) pairs")
    print(f"  Energy range: [{eigenvalues.min():.2f}, {eigenvalues.max():.2f}] eV")
    print(f"  Projectability range: [{projectabilities.min():.4f}, {projectabilities.max():.4f}]")
    scdm_mu, scdm_sigma, mu_fit, sigma_fit = fit_scdm_parameters(eigenvalues, projectabilities)
    print(f"  Fit results: mu_fit = {mu_fit:.4f} eV, sigma_fit = {sigma_fit:.4f} eV")
    print(f"  SCDM parameters: scdm_mu = {scdm_mu:.4f} eV, scdm_sigma = {scdm_sigma:.4f} eV")
    return {
        'scdm_mu': scdm_mu,
        'scdm_sigma': scdm_sigma,
        'mu_fit': mu_fit,
        'sigma_fit': sigma_fit,
    }

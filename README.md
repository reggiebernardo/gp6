# gp6

Python code for Gaussian Processes (for Physics) specifically designed for the reconstruction of late time cosmological data (e.g., [arXiv:2105.12970](https://arxiv.org/abs/2105.12970), [arXiv:2106.08688](https://arxiv.org/abs/2106.08688) and [arXiv:2111.08289](https://arxiv.org/abs/2111.08289)).

Please cite the above papers when using `gp6`, and let me know about any questions or comments. Always happy to discuss. Thanks. - **Reggie**

Installation: `pip install gp6`

##### *Minimal example* (Hubble Expansion Rate Reconstruction with Cosmic Chronometers without Covariance, in terminal, `python ex1_minimal_cc.py`): <br />

<table class="image" align="center" width="50%">
<tr><td><img src="https://github.com/reggiebernardo/gp6/blob/12045a545aec034a0887877146ec2f7defcd238f/Hz_CC_bygp6.png"></td></tr>
<tr><td class="caption">Hubble function reconstruction with CC data (No Covariance)</td></tr>
</table>

<table class="image" align="center" width="50%">
<tr><td><img src="https://github.com/reggiebernardo/gp6/blob/12045a545aec034a0887877146ec2f7defcd238f/dHdz_CC_bygp6.png"></td></tr>
<tr><td class="caption">Hubble function derivative reconstruction with CC data (No Covariance)</td></tr>
</table>

##### *Another minimal example* ($H(z)$ Reconstruction with CC Covariance + Bonus Quintessence potential and DE Equation of State GP inferrences): `quick_example.ipynb`.

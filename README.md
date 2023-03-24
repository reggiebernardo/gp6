# gp6

Python code for Gaussian Processes (for Physics) specifically designed for the reconstruction of late time cosmological data (e.g., [arXiv:2105.12970](https://arxiv.org/abs/2105.12970), [arXiv:2106.08688](https://arxiv.org/abs/2106.08688) and [arXiv:2111.08289](https://arxiv.org/abs/2111.08289)).

Please cite the above papers when using `gp6`, and let me know about any questions or comments. Always happy to discuss. Thanks. - **Reggie**

Installation: `pip install gp6`

##### *Minimal example* (Hubble Expansion Rate Reconstruction with Cosmic Chronometers without Covariance, in terminal, `python ex1_minimal_cc.py`): <br />

<table class="image" align="center" width="50%">
<tr><td><img src="./Hz_CCbygp6.png"></td></tr>
<tr><td class="caption">Hubble function reconstruction with CC data (No Covariance)</td></tr>
</table>

<table class="image" align="center" width="50%">
<tr><td><img src="./dHdz_CCbygp6.png"></td></tr>
<tr><td class="caption">Hubble function derivative reconstruction with CC data (No Covariance)</td></tr>
</table>

##### *Another minimal example* ($H(z)$ Reconstruction with CC Covariance + Bonus Quintessence potential reconstruction and DE Equation of State): `quick_example.ipynb`.

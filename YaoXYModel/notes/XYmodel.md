## XYmodel
Consider the following Hamiltonian describing a chain of N spin-1/2
$$H=\sum_{ i}J_{xy}(S_i^xS_{i+1}^x+S_i^yS_{i+1}^y)+\sum_i h S_i^z$$
where $S_i^x=\frac{1}{2}X_i,S_i^y=\frac{1}{2}Y_i,S_i^z=\frac{1}{2}Z_i$.
Denote $a^\dagger_i=S^x_i+iS^y_i$ and $a_i=S^x_i-iS^y_i$, then $S^x_i=\frac{1}{2}(a^\dagger_i+a_i)$ and $S^y_i=\frac{1}{2i}(a^\dagger_i-a_i)$.
$$a^\dagger_ia_i=(S^x_i+iS^y_i)(S^x_i-iS^y_i)=\frac{1}{2}+S^z_i$$
$$a_ia^\dagger_i=(S^x_i-iS^y_i)(S^x_i+iS^y_i)=\frac{1}{2}-S^z_i$$
$$
\begin{aligned}
H & =\sum _{i}\frac{1}{4} J_{xy}(a_{i}^{\dagger } +a_{i} )(a_{i+1}^{\dagger } +a_{i+1} )-\frac{1}{4} J_{xy} (a_{i}^{\dagger } -a_{i} )(a_{i+1}^{\dagger } -a_{i+1} )+\sum _{i}  ha_{i}^{\dagger }a_{i} +\frac{N h}{2}\\
 & = \frac{1}{2} \sum _{i}  J_{xy}(a_{i}^{\dagger } a_{i+1} +a_{i} a_{i+1}^{\dagger } )+\sum _{i}  ha_{i}^{\dagger } a_{i} +\frac{N h}{2}
\end{aligned}
$$
Jordan-Wigner transformation:
$$
c_i=(\prod^{i-1}_{j=1}Z_j)a_i, c^\dagger_i=(\prod^{i-1}_{j=1}Z_j)a^\dagger_i
$$

We have $a^\dagger_ia_i=c^\dagger_ic_i$ and
$$
c^\dagger_ic_{i+1}= a^\dagger_iZ_{i}a_{i+1}=(S^x_i+iS^y_i)Z_{i}a_{i+1}=(-iS^y_i-S^x_i)a_{i+1}=-a^\dagger_ia_{i+1}
$$
$$
c^\dagger_{i+1}c_{i}= a^\dagger_{i+1}Z_{i}a_{i}=a^\dagger_{i+1}Z_{i}(S^x_i-iS^y_i)=a^\dagger_{i+1}(iS^y_i-S^x_i)=-a^\dagger_{i+1}a_{i}
$$
Therefore $
a^\dagger_ia_{i+1}=-c^\dagger_ic_{i+1}$ and $a^\dagger_{i+1}a_i=-c^\dagger_{i+1}c_{i}$.
$$
H=-\frac{1}{2}\sum_iJ_{xy}(c^\dagger_ic_{i+1}+c^\dagger_{i+1}c_i)+\sum_i hc^\dagger_ic_i+\frac{Nh}{2}
$$
We ignore the constant term $\frac{Nh}{2}$ and denote $H$ in matrix form:
$$
H= c^\dagger A c
$$
where 
$$c=\begin{pmatrix}
c_1\\c_2\\ \vdots \\ c_N
\end{pmatrix},
A=\left(\begin{array}{ccccc}
h & -\frac{1}{2}J_x & 0 & \ldots & 0 \\
-\frac{1}{2}J_x & h & -\frac{1}{2}J_x & \ldots & 0 \\
0 & -\frac{1}{2}J_x & h & \ldots & 0 \\
\vdots & \vdots & \vdots & \ddots & \vdots \\
0 & 0 & 0 & -\frac{1}{2}J_x & h
\end{array}\right).$$
For periodic boundary condition, we have 
$$
A=\left(\begin{array}{ccccc}
h & -\frac{1}{2}J_x & 0 & \ldots & -\frac{1}{2}J_x \\
-\frac{1}{2}J_x & h & -\frac{1}{2}J_x & \ldots & 0 \\
0 & -\frac{1}{2}J_x & h & \ldots & 0 \\
\vdots & \vdots & \vdots & \ddots & \vdots \\
-\frac{1}{2}J_x & 0 & 0 & -\frac{1}{2}J_x & h
\end{array}\right) $$
Suppose the eigendecomposition of $A$ is $A=U\Lambda U^\dagger$, then we can transform $c$ into the momentum space $c=U\eta$ and the Hamiltonian becomes
$$
H=c^\dagger A c=\eta^\dagger U^\dagger U\Lambda U^\dagger U\eta=\eta^\dagger \Lambda \eta = \sum_k\Lambda_k\eta^\dagger_k\eta_k
$$
Also, we have $\eta=U^\dagger c$ and
$$
\eta_k=\sum_jU^*_{jk}c_j, \eta^\dagger_k=\sum_jU_{jk}c^\dagger_j
$$

## Floquet Hamiltonian
We apply a periodic magnetic field in the $x$ direction on the first spin, the Hamiltonian of this magnetic field is 
$$H_1=\Omega \cos(\omega t)S^x_1 = \frac{1}{2}\Omega (\cos(\omega t)S^x_1+\sin(\omega t)S^y_1)+\frac{1}{2}\Omega (\cos(\omega t)S^x_1-\sin(\omega t)S^y_1)$$ 

According to the rotating-wave approximation, we can omit the second term. The Hamiltonian becomes
$$
\begin{aligned}H_1&= \frac{1}{2}\Omega (\cos(\omega t)S^x_1+\sin(\omega t)S^y_1)\\
&=\frac{1}{4}\Omega (\cos(\omega t)(c^\dagger_1+c_1)+\frac{1}{i}\sin(\omega t)(c^\dagger_1-c_1))\\
&=\frac{1}{4}\Omega \cos(\omega t)\sum_{k}(U_{1k}^*\eta^\dagger_k+U_{1k}\eta_k)+\frac{1}{4i}\Omega \sin(\omega t)\sum_{k}(U_{1k}^*\eta^\dagger_k-U_{1k}\eta_k)\\
&=\frac{1}{4}\Omega \cos(\omega t)\sum_{k}((r_{k}-is_k)\eta^\dagger_k+(r_{k}+is_k)\eta_k)+\frac{1}{4i}\Omega \sin(\omega t)\sum_{k}((r_{k}-is_k)\eta^\dagger_k-(r_{k}+is_k)\eta_k)\\
&=\frac{1}{4}\Omega \sum_{k}(r_{k}\cos(\omega t)-s_k\sin(\omega t))(\eta^\dagger_k+\eta_k)+\frac{1}{4i}\Omega \sum_{k}(r_{k}\sin(\omega t)+s_k\cos(\omega t))(\eta^\dagger_k-\eta_k)\\
\end{aligned}$$ 
where $r_k$ and $s_k$ are the real and imaginary parts of $U_{1k}$.
We focus on a single mode $k$, regard the subsystem as a two-level system. Then $(\eta^\dagger_k+\eta_k)$ is the $X$ operator and $i(\eta^\dagger_k-\eta_k)$ is the $Y$ operator. The Hamiltonian becomes
$$
\begin{aligned}
H_{k1}&=\frac{1}{4}\Omega (r_{k}\cos(\omega t)-s_k\sin(\omega t))X-\frac{1}{4}\Omega (r_{k}\sin(\omega t)+s_k\cos(\omega t))Y\\
&=\frac{1}{4}\Omega \sqrt{r^2_k+s^2_k}\cos(\omega t+\phi_k)X-\frac{1}{4}\Omega \sqrt{r^2_k+s^2_k}\sin(\omega t+\phi_k)Y
\end{aligned}
$$
where $\cos(\phi_k)=r_k/\sqrt{r^2_k+s^2_k}$ and $\sin(\phi_k)=s_k/\sqrt{r^2_k+s^2_k}$.
And $
H_{k}=\Lambda_k \eta_k^\dagger \eta_k=\frac{1}{2}\Lambda_kZ
$ (ingore the constant term $\frac{1}{2}\Lambda_k$).
Toy model of two subsystems with an anti-commute term in the Hamiltonian:
$$
H=\Delta_1 Z_1+\Omega_1X_1+\Delta_2 Z_2+\Omega_2Z_1X_2
$$

## Rotation frame(Nilsen and Chuang, 7.7.2)
We can denote the two level system hamiltonian $H_k+H_{k1}$ by
$$
H_{\text{two level}}= \frac{\omega_0}{2}Z+g(X\cos(\omega t)-Y\sin(\omega t))
$$
where 
$$\omega_0 =\Lambda_k ,g=\frac{1}{4}\Omega \sqrt{r^2_k+s^2_k}.$$
The state of the system is $|\chi(t)\rangle$ and define the state in the rotation frame as $|\varphi(t)\rangle=e^{-i \omega t Z / 2}|\chi(t)\rangle$, such that the Schrödinger equation
$$
i \partial_t|\chi(t)\rangle=H_{\text{two level}}|\chi(t)\rangle
$$
can be re-expressed as
$$
i \partial_t|\varphi(t)\rangle=\left[e^{-i \omega Z t / 2} H_{\text{two level}} e^{i \omega Z t / 2}+\frac{\omega}{2} Z\right]|\varphi(t)\rangle .
$$

Since
$$
e^{-i \omega Z t / 2} X e^{i \omega Z t / 2}=(X \cos \omega t+Y \sin \omega t)
$$
and 
$$
e^{-i \omega Z t / 2} Y e^{i \omega Z t / 2}=(-X \sin \omega t+Y \cos \omega t)
$$
we have
$$
i \partial_t|\varphi(t)\rangle=\left[\frac{\omega_0+\omega}{2} Z+g X\right]|\varphi(t)\rangle,
$$
where the terms on the right multiplying the state can be identified as the effective rotation frame Hamiltonian. The solution to this equation is
$$
|\varphi(t)\rangle=e^{-i\left[\frac{\omega_0+\omega}{2} Z+g X\right] t}|\varphi(0)\rangle .
$$

The concept of resonance arises from the behavior of this solution, which can be understood to be a single qubit rotation about the axis
$$
\hat{n}=\frac{\hat{z}+\frac{2 g}{\omega_0+\omega} \hat{x}}{\sqrt{1+\left(\frac{2 g}{\omega_0+\omega}\right)^2}}
$$
by an angle
$$
|\vec{n}|=t \sqrt{\left(\frac{\omega_0+\omega}{2}\right)^2+g^2}.
$$
Therefore,
$$|\chi(t)\rangle=e^{i \omega t Z / 2}|\varphi(t)\rangle =e^{i \omega t Z / 2}e^{-i\left[\frac{\omega_0+\omega}{2} Z+g X\right] t}|\varphi(0)\rangle=e^{i \omega t Z / 2}e^{-i\left[\frac{\omega_0+\omega}{2} Z+g X\right] t}|\chi(0)\rangle$$
## Energy level with periodic boundary condition
For the periodic boundary condition, the transformation matrix is the quantum fourier transformation:
$$
\eta_k=\frac{1}{\sqrt{N}}\sum_{j}e^{2\pi ikj/N}c_j,
\eta^\dagger_k=\frac{1}{\sqrt{N}}\sum_{j}e^{-2\pi ikj/N}c^\dagger_j.
$$
$$
c_j=\frac{1}{\sqrt{N}}\sum_{k}e^{-2\pi ikj/N}\eta_k,
c^\dagger_j=\frac{1}{\sqrt{N}}\sum_{k}e^{2\pi ikj/N}\eta^\dagger_k.
$$
$$
c^\dagger_jc_{j+1}=\frac{1}{N}\sum_{k_1,k_2}e^{2\pi i(k_1j-k_2j-k_2)/N}\eta^\dagger_{k_1}\eta_{k_2},
$$
$$
c^\dagger_{j+1}c_j=\frac{1}{N}\sum_{k_1,k_2}e^{2\pi i(k_1j-k_2j+k_1)/N}\eta^\dagger_{k_1}\eta_{k_2},
$$
$$
c^\dagger_jc_j=\frac{1}{N}\sum_{k_1,k_2}e^{2\pi i(k_1-k_2)j/N}\eta^\dagger_{k_1}\eta_{k_2},
$$
$$
\begin{aligned}
H&=\frac{1}{2}\sum_j(c^\dagger_jc_{j+1}+c^\dagger_{j+1}c_j)+\sum_j c^\dagger_jc_j+\frac{N}{2}\\
&= \frac{1}{2}\sum_j(\frac{1}{N}\sum_{k_1,k_2}e^{2\pi i(k_1j-k_2j-k_2)/N}\eta^\dagger_{k_1}\eta_{k_2}+\frac{1}{N}\sum_{k_1,k_2}e^{2\pi i(k_1j-k_2j+k_1)/N}\eta^\dagger_{k_1}\eta_{k_2})+\sum_j \frac{1}{N}\sum_{k_1,k_2}e^{2\pi i(k_1-k_2)j/N}\eta^\dagger_{k_1}\eta_{k_2}+\frac{N}{2}\\
&=\frac{1}{2N}\sum_{j,k_1,k_2}e^{2\pi i(k_1-k_2)j/N}(\eta^\dagger_{k_1}\eta_{k_2}e^{-2\pi ik_2/N}+\eta^\dagger_{k_1}\eta_{k_2}e^{2\pi ik_1/N})+\frac{1}{N}\sum_{j,k_1,k_2}e^{2\pi i(k_1-k_2)j/N}\eta^\dagger_{k_1}\eta_{k_2}+\frac{N}{2}\\
&=\frac{1}{2N}\sum_{k_1,k_2}(\eta^\dagger_{k_1}\eta_{k_2}e^{-2\pi ik_2/N}+\eta^\dagger_{k_1}\eta_{k_2}e^{2\pi ik_1/N})\sum_{j}e^{2\pi i(k_1-k_2)j/N}+\frac{1}{N}\sum_{k_1,k_2}\eta^\dagger_{k_1}\eta_{k_2}\sum_{j}e^{2\pi i(k_1-k_2)j/N}+\frac{N}{2}\\
&=\frac{1}{2N}\sum_{k_1,k_2}(\eta^\dagger_{k_1}\eta_{k_2}e^{-2\pi ik_2/N}+\eta^\dagger_{k_1}\eta_{k_2}e^{2\pi ik_1/N})N\delta_{k_1,k_2}+\frac{1}{N}\sum_{k_1,k_2}\eta^\dagger_{k_1}\eta_{k_2}N\delta_{k_1,k_2}+\frac{N}{2}\\
&=\frac{1}{2}\sum_{k}(\eta^\dagger_{k}\eta_{k}e^{-2\pi ik/N}+\eta^\dagger_{k}\eta_{k}e^{2\pi ik/N})+\sum_{k}\eta^\dagger_{k}\eta_{k}+\frac{N}{2}\\
&=\sum_{k}(\cos(2πk/N)+1)\eta^\dagger_{k}\eta_{k}+\frac{N}{2}
\end{aligned}
$$
The ground state of $H$ is the vacuum state of $\eta_k$ and $ | 11...1 \rangle $ in the original spin language. The ground state energy is $E_0=N/2$.

syms x y
N0_L0=3;
m_H=4;
mu0_L0=3;
c0_c1=2
a_H=2
b_H=4
P_mu=2
a11=i*sqrt(y-N0_L0)-m_H/x;
a12=-(i*sqrt(y-N0_L0)+m_H/x)

z=(-2*x)*(1/a_H)*(b_H/a_H)*y*(c0_c1)+2*x*(1/a_H);
R=((1-(b_H/a_H)*y*c0_c1)^(-1/2))*(0.25*P_mu*(x/a_H)+0.5*y*(c0_c1)*(x/a_H)-0.5*(b_H/a_H)*(x/a_H)*y*(c0_c1));
WW=(exp(-z/2))*(z^R)*(1-(1/z)*(0.5-R)^2+((0.5-R)^2)*((1.5-R)^2)*(1/2)*(1/z^2))
DWW=z^R*exp(-z/2)*((R - 1/2)^2/z^2 - ((R - 1/2)^2*(R - 3/2)^2)/z^3) - (z^R*exp(-z/2)*(((R - 1/2)^2*(R - 3/2)^2)/(2*z^2) - (R - 1/2)^2/z + 1))/2 + R*z^(R - 1)*exp(-z/2)*(((R - 1/2)^2*(R - 3/2)^2)/(2*z^2) - (R - 1/2)^2/z + 1);
a13=-(mu0_L0)*((1/2)*(a_H/x)*WW+DWW)

a21=(i*sqrt(y-N0_L0)-i*m_H*sqrt(y-N0_L0)-m_H/x)*exp(-i*x*sqrt(y-N0_L0));
a22=-(i*sqrt(y-N0_L0)-i*m_H*sqrt(y-N0_L0)+m_H/x)*exp(i*x*sqrt(y-N0_L0));
a23=0;

a31=1;
a32=1;
a33=WW;

AA=[a11 a12 a13; a21 a22 a23; a31 a32 a33]

exp1=real(det(AA))


hold on
N0_L0=3;
m_H=8;
mu0_L0=3;
c0_c1=2
a_H=2
b_H=4
P_mu=2
a11=i*sqrt(y-N0_L0)-m_H/x;
a12=-(i*sqrt(y-N0_L0)+m_H/x)

z=(-2*x)*(1/a_H)*(b_H/a_H)*y*(c0_c1)+2*x*(1/a_H);
R=((1-(b_H/a_H)*y*c0_c1)^(-1/2))*(0.25*P_mu*(x/a_H)+0.5*y*(c0_c1)*(x/a_H)-0.5*(b_H/a_H)*(x/a_H)*y*(c0_c1));
WW=(exp(-z/2))*(z^R)*(1-(1/z)*(0.5-R)^2+((0.5-R)^2)*((1.5-R)^2)*(1/2)*(1/z^2))  %WhittakerW expansion at z=0
DWW=z^R*exp(-z/2)*((R - 1/2)^2/z^2 - ((R - 1/2)^2*(R - 3/2)^2)/z^3) - (z^R*exp(-z/2)*(((R - 1/2)^2*(R - 3/2)^2)/(2*z^2) - (R - 1/2)^2/z + 1))/2 + R*z^(R - 1)*exp(-z/2)*(((R - 1/2)^2*(R - 3/2)^2)/(2*z^2) - (R - 1/2)^2/z + 1);  %Differentiation of whittakerW at z=0
a13=-(mu0_L0)*((1/2)*(a_H/x)*WW+DWW)

a21=(i*sqrt(y-N0_L0)-i*m_H*sqrt(y-N0_L0)-m_H/x)*exp(-i*x*sqrt(y-N0_L0));
a22=-(i*sqrt(y-N0_L0)-i*m_H*sqrt(y-N0_L0)+m_H/x)*exp(i*x*sqrt(y-N0_L0));
a23=0;

a31=1;
a32=1;
a33=WW;

AA=[a11 a12 a13; a21 a22 a23; a31 a32 a33]
exp2=real(det(AA))


axis([1.5,6,1,2.8]);


ezplot(exp1);
hold on
ezplot(exp2);


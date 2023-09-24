(* ::Package:: *)

BeginPackage["OneDfunctional`"]

betafunctional::usage = "betafunctional[\!\(\*SubscriptBox[\(\[CapitalDelta]\), \(\[Phi]\)]\), m, \[CapitalDelta]] returns \!\(\*SubscriptBox[\(\[Beta]\), \(m\)]\)(\[CapitalDelta]) at external dimension \!\(\*SubscriptBox[\(\[CapitalDelta]\), \(\[Phi]\)]\).";
betafunctionalD2::usage = "betafunctional[\!\(\*SubscriptBox[\(\[CapitalDelta]\), \(\[Phi]\)]\), m, \[CapitalDelta]] returns \!\(\*SubscriptBox[\(d2\[Beta]\), \(m\)]\)(\[CapitalDelta]) when \[CapitalDelta]=2\[CapitalDelta]\[Phi]+2m at external dimension \!\(\*SubscriptBox[\(\[CapitalDelta]\), \(\[Phi]\)]\).";

betafunctionallist::usage = "betafunctionallist[\!\(\*SubscriptBox[\(\[CapitalDelta]\), \(\[Phi]\)]\), m, \[CapitalDelta]] returns {\!\(\*SubscriptBox[\(\[Beta]\), \(i\)]\)(\[CapitalDelta]) (i=1,2,..m)} at external dimension \!\(\*SubscriptBox[\(\[CapitalDelta]\), \(\[Phi]\)]\).";
alphafunctional::usage = "alphafunctional[\!\(\*SubscriptBox[\(\[CapitalDelta]\), \(\[Phi]\)]\), m, \[CapitalDelta]] returns \!\(\*SubscriptBox[\(\[Alpha]\), \(m\)]\)(\[CapitalDelta]) at external dimension \!\(\*SubscriptBox[\(\[CapitalDelta]\), \(\[Phi]\)]\).";
alphafunctionallist::usage = "alphafunctionallist[\!\(\*SubscriptBox[\(\[CapitalDelta]\), \(\[Phi]\)]\), m, \[CapitalDelta]] returns {\!\(\*SubscriptBox[\(\[Alpha]\), \(i\)]\)(\[CapitalDelta]) (i=0,1,..m)} at external dimension \!\(\*SubscriptBox[\(\[CapitalDelta]\), \(\[Phi]\)]\).";
functional::usage = "functional[\!\(\*SubscriptBox[\(\[CapitalDelta]\), \(\[Phi]\)]\), m, \[CapitalDelta]] returns {\!\(\*SubscriptBox[\(\[Beta]\), \(m\)]\)(\[CapitalDelta]), \!\(\*SubscriptBox[\(\[Alpha]\), \(m\)]\)(\[CapitalDelta])} at external dimension \!\(\*SubscriptBox[\(\[CapitalDelta]\), \(\[Phi]\)]\).";
functionallist::usage = "functional[\!\(\*SubscriptBox[\(\[CapitalDelta]\), \(\[Phi]\)]\), m, \[CapitalDelta]] returns {betafunctionallist[\!\(\*SubscriptBox[\(\[CapitalDelta]\), \(\[Phi]\)]\), m, \[CapitalDelta]],  alphafunctionallist[\!\(\*SubscriptBox[\(\[CapitalDelta]\), \(\[Phi]\)]\), m, \[CapitalDelta]]}.";
betatildezero::usage = "betatildezero[\!\(\*SubscriptBox[\(\[CapitalDelta]\), \(\[Phi]\)]\), \[CapitalDelta]] returns \!\(\*SubscriptBox[OverscriptBox[\(\[Beta]\), \(~\)], \(0\)]\)[\!\(\*SubscriptBox[\(\[CapitalDelta]\), \(\[Phi]\)]\), \[CapitalDelta]]."
tilde0D2N::usage = "tilde0D2N[a,n] returns \!\(\*SuperscriptBox[SubscriptBox[OverscriptBox[\(\[Beta]\), \(~\)], \(0\)], \(\*\"\\\"\<\>\"\)]\) at extdim a and dimension 2a+2n for n strict positive integer";
betafunctionalderi2::usage="betafunctionalderi2[a,m,n] returns \!\(\*SuperscriptBox[SubscriptBox[\(\[Beta]\), \(m\)], \(\*\"\\\"\<\>\"\)]\) at extdim a and dimension 2a+2n for n strict positive integer";
alphatildefunctional::usage="alphatildefunctional[a,m,\[CapitalDelta]] returns \!\(\*SubscriptBox[OverscriptBox[\(\[Alpha]\), \(~\)], \(m\)]\)(\[CapitalDelta]).";
betatildefunctional::usage="betatildefunctional[a,m,\[CapitalDelta]] returns \!\(\*SubscriptBox[OverscriptBox[\(\[Beta]\), \(~\)], \(m\)]\)(\[CapitalDelta]).";
cn::usage="cn[\!\(\*SubscriptBox[\(\[CapitalDelta]\), \(\[Phi]\)]\)][n] return \!\(\*SubscriptBox[\(c\), \(n\)]\).";
dn::usage="dn[\!\(\*SubscriptBox[\(\[CapitalDelta]\), \(\[Phi]\)]\)][n] return \!\(\*SubscriptBox[\(d\), \(n\)]\).";
afermion::usage="afermion[\!\(\*SubscriptBox[\(\[CapitalDelta]\), \(\[Phi]\)]\)][m] return generalized free fermion OPE^2 \!\(\*SubscriptBox[\(a\), \(m\)]\) at 2\!\(\*SubscriptBox[\(\[CapitalDelta]\), \(\[Phi]\)]\)+2m+1.";
am::usage="am[\!\(\*SubscriptBox[\(\[CapitalDelta]\), \(\[Phi]\)]\)][m] return generalized free boson OPE^2 \!\(\*SubscriptBox[\(a\), \(m\)]\) at 2\!\(\*SubscriptBox[\(\[CapitalDelta]\), \(\[Phi]\)]\)+2m.";
cb1::usage="cb1[\!\(\*SubscriptBox[\(\[CapitalDelta]\), \(\[Phi]\)]\)][\[CapitalDelta],z] return 1d Conformal block.";
functionallistPara::usage="";


Begin["`Private`"]



bs0[n_]=(Sqrt[Pi]*Gamma[n+\[CapitalDelta]\[Phi]]^4*Gamma[-1/2+n+2*\[CapitalDelta]\[Phi]]^2)/((2*n-\[CapitalDelta]+2*\[CapitalDelta]\[Phi])*(-1+2*n+\[CapitalDelta]+2*\[CapitalDelta]\[Phi])*n!^2*Gamma[\[CapitalDelta]\[Phi]]^4*Gamma[2*(n+\[CapitalDelta]\[Phi])]*Gamma[-1/2+2*n+2*\[CapitalDelta]\[Phi]]);
as0[n_]=((-1)^(1)*Sqrt[Pi]*Gamma[n+\[CapitalDelta]\[Phi]]^4*Gamma[-1/2+n+2*\[CapitalDelta]\[Phi]]^2*(-1+4*n+4*\[CapitalDelta]\[Phi]+(2*n-\[CapitalDelta]+2*\[CapitalDelta]\[Phi])*(-1+2*n+\[CapitalDelta]+2*\[CapitalDelta]\[Phi])*(HarmonicNumber[n]-2*HarmonicNumber[-1+n+\[CapitalDelta]\[Phi]]-HarmonicNumber[-3/2+n+2*\[CapitalDelta]\[Phi]]+2*HarmonicNumber[-2+4*n+4*\[CapitalDelta]\[Phi]]-Log[4])))/((-2*n+\[CapitalDelta]-2*\[CapitalDelta]\[Phi])^2*(-1+2*n+\[CapitalDelta]+2*\[CapitalDelta]\[Phi])^2*n!^2*Gamma[\[CapitalDelta]\[Phi]]^4*Gamma[2*(n+\[CapitalDelta]\[Phi])]*Gamma[-1/2+2*n+2*\[CapitalDelta]\[Phi]]);
A[\[CapitalDelta]_,\[CapitalDelta]\[Phi]_]=(4^(-2+\[CapitalDelta])*Gamma[\[CapitalDelta]/2]^4*Gamma[(-1+\[CapitalDelta])/2+\[CapitalDelta]\[Phi]]^2*Gamma[-1/2*\[CapitalDelta]+\[CapitalDelta]\[Phi]]^2)/(Gamma[2*\[CapitalDelta]]*Gamma[\[CapitalDelta]\[Phi]]^4);
AD2[\[CapitalDelta]_,\[CapitalDelta]\[Phi]_]=(4^(-2+\[CapitalDelta])*Gamma[\[CapitalDelta]/2]^4*Gamma[(-1+\[CapitalDelta])/2+\[CapitalDelta]\[Phi]]^2)/(Gamma[2*\[CapitalDelta]]*Gamma[\[CapitalDelta]\[Phi]]^4);

bs0OverA[n_,\[CapitalDelta]_,\[CapitalDelta]\[Phi]_]=(2^(3-2 \[CapitalDelta]) Sqrt[\[Pi]] Gamma[2 \[CapitalDelta]] Gamma[1+\[CapitalDelta]/2-\[CapitalDelta]\[Phi]]^2 Gamma[n+\[CapitalDelta]\[Phi]]^4 Gamma[-(1/2)+n+2 \[CapitalDelta]\[Phi]]^2)/((-1+2 n+\[CapitalDelta]+2 \[CapitalDelta]\[Phi]) (n!)^2 Gamma[\[CapitalDelta]/2]^4 Gamma[1-n+\[CapitalDelta]/2-\[CapitalDelta]\[Phi]]^2 Gamma[2 (n+\[CapitalDelta]\[Phi])] Gamma[1/2 (-1+\[CapitalDelta])+\[CapitalDelta]\[Phi]]^2 Gamma[n-\[CapitalDelta]/2+\[CapitalDelta]\[Phi]] Gamma[1+n-\[CapitalDelta]/2+\[CapitalDelta]\[Phi]] Gamma[-(1/2)+2 n+2 \[CapitalDelta]\[Phi]]);
as0OverA[n_,\[CapitalDelta]_,\[CapitalDelta]\[Phi]_]=-((4^(1-\[CapitalDelta]) Sqrt[\[Pi]] Gamma[2 \[CapitalDelta]] Gamma[1+\[CapitalDelta]/2-\[CapitalDelta]\[Phi]]^2 Gamma[n+\[CapitalDelta]\[Phi]]^4 Gamma[-(1/2)+n+2 \[CapitalDelta]\[Phi]]^2 (-1+4 n+4 \[CapitalDelta]\[Phi]+(2 n-\[CapitalDelta]+2 \[CapitalDelta]\[Phi]) (-1+2 n+\[CapitalDelta]+2 \[CapitalDelta]\[Phi]) (HarmonicNumber[n]-2 HarmonicNumber[-1+n+\[CapitalDelta]\[Phi]]-HarmonicNumber[-(3/2)+n+2 \[CapitalDelta]\[Phi]]+2 HarmonicNumber[-2+4 n+4 \[CapitalDelta]\[Phi]]-Log[4])))/((-1+2 n+\[CapitalDelta]+2 \[CapitalDelta]\[Phi])^2 (n!)^2 Gamma[\[CapitalDelta]/2]^4 Gamma[1-n+\[CapitalDelta]/2-\[CapitalDelta]\[Phi]]^2 Gamma[2 (n+\[CapitalDelta]\[Phi])] Gamma[1/2 (-1+\[CapitalDelta])+\[CapitalDelta]\[Phi]]^2 Gamma[1+n-\[CapitalDelta]/2+\[CapitalDelta]\[Phi]]^2 Gamma[-(1/2)+2 n+2 \[CapitalDelta]\[Phi]]));


bCS[\[CapitalDelta]\[Phi]_,n_]=(2^(-2 (n+\[CapitalDelta]+\[CapitalDelta]\[Phi])) Sqrt[\[Pi]] Gamma[2 \[CapitalDelta]] Gamma[n+\[CapitalDelta]\[Phi]]^3 ((Sqrt[\[Pi]] (\[CapitalDelta]-2 \[CapitalDelta]\[Phi]) (-1+\[CapitalDelta]+2 \[CapitalDelta]\[Phi]))/(Gamma[1-\[CapitalDelta]/2+\[CapitalDelta]\[Phi]]^2 Gamma[(1+\[CapitalDelta])/2+\[CapitalDelta]\[Phi]]^2)-(32 \[Beta]0s0 Gamma[2 \[CapitalDelta]\[Phi]])/(Gamma[1/2 (-1+\[CapitalDelta])+\[CapitalDelta]\[Phi]]^2 Gamma[-(\[CapitalDelta]/2)+\[CapitalDelta]\[Phi]]^2 Gamma[-(1/2)+2 \[CapitalDelta]\[Phi]])) Gamma[-(1/2)+n+2 \[CapitalDelta]\[Phi]]^2)/(Gamma[1+n]^2 Gamma[\[CapitalDelta]/2]^4 Gamma[1/2+n+\[CapitalDelta]\[Phi]] Gamma[-(1/2)+2 n+2 \[CapitalDelta]\[Phi]]);
bCSD2[\[CapitalDelta]\[Phi]_,n_]=((-1+4*\[CapitalDelta]\[Phi])*Gamma[2*\[CapitalDelta]]*Gamma[\[CapitalDelta]\[Phi]]*Gamma[1/2+\[CapitalDelta]\[Phi]]*Gamma[n+\[CapitalDelta]\[Phi]]^3*(2*\[Beta]0s0*(\[CapitalDelta]-\[CapitalDelta]^2+2*\[CapitalDelta]\[Phi]*(-1+2*\[CapitalDelta]\[Phi]))*Gamma[2*\[CapitalDelta]\[Phi]]+Sqrt[Pi]*Gamma[-1/2+2*\[CapitalDelta]\[Phi]])*Gamma[-1/2+n+2*\[CapitalDelta]\[Phi]]^2)/(2^(2*(-1+n+\[CapitalDelta]))*(\[CapitalDelta]-2*\[CapitalDelta]\[Phi])*(-1+\[CapitalDelta]+2*\[CapitalDelta]\[Phi])*Gamma[1+n]^2*Gamma[\[CapitalDelta]/2]^4*Gamma[2*\[CapitalDelta]\[Phi]]*Gamma[1/2+n+\[CapitalDelta]\[Phi]]*Gamma[(-1+\[CapitalDelta])/2+\[CapitalDelta]\[Phi]]^2*Gamma[1/2+2*\[CapitalDelta]\[Phi]]*Gamma[-1/2+2*n+2*\[CapitalDelta]\[Phi]]);

aCS[\[CapitalDelta]\[Phi]_,n_]=(2^(1-2 n-2 \[CapitalDelta]) (-1+4 \[CapitalDelta]\[Phi]) Gamma[2 \[CapitalDelta]] Gamma[\[CapitalDelta]\[Phi]] Gamma[1/2+\[CapitalDelta]\[Phi]] Gamma[n+\[CapitalDelta]\[Phi]]^3 (2 \[Beta]0s0 (\[CapitalDelta]-\[CapitalDelta]^2+2 \[CapitalDelta]\[Phi] (-1+2 \[CapitalDelta]\[Phi])) Gamma[2 \[CapitalDelta]\[Phi]]+Sqrt[\[Pi]] Gamma[-(1/2)+2 \[CapitalDelta]\[Phi]]) Gamma[-(1/2)+n+2 \[CapitalDelta]\[Phi]]^2 (Gamma[1+\[CapitalDelta]\[Phi]] Log[4]+\[CapitalDelta]\[Phi] Gamma[\[CapitalDelta]\[Phi]] (2 PolyGamma[0,1+n]-3 PolyGamma[0,n+\[CapitalDelta]\[Phi]]+PolyGamma[0,1/2+n+\[CapitalDelta]\[Phi]]-2 PolyGamma[0,-(1/2)+n+2 \[CapitalDelta]\[Phi]]+2 PolyGamma[0,-(1/2)+2 n+2 \[CapitalDelta]\[Phi]])))/((-1+\[CapitalDelta]+2 \[CapitalDelta]\[Phi]) Gamma[1+n]^2 Gamma[\[CapitalDelta]/2]^4 Gamma[2 \[CapitalDelta]\[Phi]] Gamma[1+\[CapitalDelta]\[Phi]] Gamma[1/2+n+\[CapitalDelta]\[Phi]] Gamma[1/2 (-1+\[CapitalDelta])+\[CapitalDelta]\[Phi]]^2 Gamma[-(\[CapitalDelta]/2)+\[CapitalDelta]\[Phi]] Gamma[1-\[CapitalDelta]/2+\[CapitalDelta]\[Phi]] Gamma[1/2+2 \[CapitalDelta]\[Phi]] Gamma[-(1/2)+2 n+2 \[CapitalDelta]\[Phi]]);
\[Beta]SS[n_]= bs0OverA[n,\[CapitalDelta],\[CapitalDelta]\[Phi]]+2* b0[2n]/A[\[CapitalDelta],\[CapitalDelta]\[Phi]]+2bCS[\[CapitalDelta]\[Phi],n]; (*first term is simplified version of bs0[n]/A[\[CapitalDelta],\[CapitalDelta]\[Phi]]*)

\[Beta]SSD2[n_]= bs0[n]/AD2[\[CapitalDelta],\[CapitalDelta]\[Phi]]+2* b0[2n]/AD2[\[CapitalDelta],\[CapitalDelta]\[Phi]]+2bCSD2[\[CapitalDelta]\[Phi],n];

\[Alpha]SS[n_]=(2*a0[2*n])/A[\[CapitalDelta],\[CapitalDelta]\[Phi]]+aCS[\[CapitalDelta]\[Phi],n]+as0OverA[n,\[CapitalDelta],\[CapitalDelta]\[Phi]];
\[Mu][n_]:=-n^2;
\[Nu][n_]:=(-1+\[CapitalDelta]) \[CapitalDelta]+\[CapitalDelta]\[Phi]+1/2 n (-1+n+4 \[CapitalDelta]\[Phi]);
\[Rho][n_]:=-(((n+2 \[CapitalDelta]\[Phi])^2 (-1+n+4 \[CapitalDelta]\[Phi])^2)/(4 (-1+2 n+4 \[CapitalDelta]\[Phi]) (1+2 n+4 \[CapitalDelta]\[Phi])));
\[Mu]p[n_]:=-2 n;
\[Nu]p[n_]:=1/2 (-1+n+2 \[CapitalDelta]\[Phi])+1/2 (n+2 \[CapitalDelta]\[Phi]);
\[Rho]p[n_]:=-(((n+2 \[CapitalDelta]\[Phi]) (-1+n+4 \[CapitalDelta]\[Phi]) (1+4 n^3-6 \[CapitalDelta]\[Phi]+24 n^2 \[CapitalDelta]\[Phi]+32 \[CapitalDelta]\[Phi]^3+n (-2+48 \[CapitalDelta]\[Phi]^2)))/(2 (-1+2 n+4 \[CapitalDelta]\[Phi])^2 (1+2 n+4 \[CapitalDelta]\[Phi])^2));
R[n_]:=(\[Pi]^(1/2) (-1)^(1-2n) Gamma[n+\[CapitalDelta]\[Phi]]^4 Gamma[-1/2+n+2\[CapitalDelta]\[Phi]]^2)/((n!)^2 Gamma[\[CapitalDelta]\[Phi]]^4 Gamma[2(n+\[CapitalDelta]\[Phi])]Gamma[2(n+\[CapitalDelta]\[Phi])-1/2]);
S[n_]:=(\[Pi]^(1/2) (-1)^(-2n) Gamma[n+\[CapitalDelta]\[Phi]]^4 Gamma[-1/2+n+2\[CapitalDelta]\[Phi]]^2)/((n!)^2 Gamma[\[CapitalDelta]\[Phi]]^4 Gamma[2(n+\[CapitalDelta]\[Phi])]Gamma[2(n+\[CapitalDelta]\[Phi])-1/2]) (-2HarmonicNumber[n+\[CapitalDelta]\[Phi]-1]-HarmonicNumber[n+2\[CapitalDelta]\[Phi]-3/2]+2HarmonicNumber[4n+4\[CapitalDelta]\[Phi]-2]+HarmonicNumber[n]-Log[4]);

defaultPrec=100;
Options[betatildezero]={PrecisionGoal->defaultPrec,sym->False};
Options[functionallist]=Options[functional]=Options[alphafunctionallist]=Options[alphafunctional]=Options[betafunctionallist]=Options[betafunctional]=Options[tilde0D2N]=Options[betafunctionalderi2]=Options[betatildefunctional]=Options[alphatildefunctional]={PrecisionGoal->defaultPrec};

betatildezero[a_,x_,OptionsPattern[]]:=If[OptionValue[sym],2^(-5-2*x)*Sqrt[Pi]*Gamma[-1+a]^3*Gamma[-5/2+2*a]*((-4096*Gamma[2*(-2+x)]*((4^(2-a)*Sqrt[Pi]*(2+(-10*a+4*a^2-(-5+x)*x)^(-1)))/(Gamma[-1/2+a]*Gamma[-5/2+a+x/2]^2)-((9+2*a-10*x+2*x^2)*Gamma[-1+a]*Gamma[-11/4+a+x/2]*Gamma[(-2+x)/2]^2*Gamma[-9/2+2*a+x]*HypergeometricPFQRegularized[{(-2+x)/2,(-2+x)/2,-3+2*a,(-5+2*a+x)/2,(-5+2*a+x)/2,-7/4+a+x/2,-11/2+2*a+x},{-3/2+x,-11/4+a+x/2,-2+a+x/2,-2+a+x/2,(-7+4*a+x)/2,(-7+4*a+x)/2},1])/2))/(Gamma[a-x/2]^2*Gamma[(-2+x)/2]^4)+(128*(-1-x+x^2)*Gamma[2*x]*((4^(2-a)*Sqrt[Pi]*(2+(6-10*a+4*a^2+x-x^2)^(-1)))/(Gamma[-1/2+a]*Gamma[-3/2+a+x/2]^2)-((-3+2*a-2*x+2*x^2)*Gamma[-1+a]*Gamma[-7/4+a+x/2]*Gamma[x/2]^2*Gamma[-5/2+2*a+x]*HypergeometricPFQRegularized[{x/2,x/2,-3+2*a,(-3+2*a+x)/2,(-3+2*a+x)/2,-3/4+a+x/2,-7/2+2*a+x},{1/2+x,-7/4+a+x/2,-1+a+x/2,-1+a+x/2,(-5+4*a+x)/2,(-5+4*a+x)/2},1])/2))/((-3+2*x)*(1+2*x)*Gamma[-1+a-x/2]^2*Gamma[x/2]^4)-(x^2*(1+x)^2*Gamma[2*(2+x)]*((4^(2-a)*Sqrt[Pi]*(2+(4-10*a+4*a^2-3*x-x^2)^(-1)))/(Gamma[-1/2+a]*Gamma[-1/2+a+x/2]^2)-((1+2*a+6*x+2*x^2)*Gamma[-1+a]*Gamma[-3/4+a+x/2]*Gamma[(2+x)/2]^2*Gamma[-1/2+2*a+x]*HypergeometricPFQRegularized[{(2+x)/2,(2+x)/2,-3+2*a,(-1+2*a+x)/2,(-1+2*a+x)/2,1/4+a+x/2,-3/2+2*a+x},{5/2+x,-3/4+a+x/2,a+x/2,a+x/2,(-3+4*a+x)/2,(-3+4*a+x)/2},1])/2))/((1+2*x)^2*(-3+4*x+4*x^2)*Gamma[-2+a-x/2]^2*Gamma[(2+x)/2]^4)),Block[{\[CapitalDelta]=SetPrecision[x,OptionValue[PrecisionGoal]],\[CapitalDelta]\[Phi]=SetPrecision[a,OptionValue[PrecisionGoal]]},betatildezero[\[CapitalDelta]\[Phi],\[CapitalDelta],sym->True]]];

betafunctional[ext_, order_, x_, OptionsPattern[]]:= Block[{
\[CapitalDelta]\[Phi] = SetPrecision[ext, OptionValue[PrecisionGoal]],
 \[CapitalDelta] = SetPrecision[x, OptionValue[PrecisionGoal]], \[Beta]0s0, b0, brule}, 
\[Beta]0s0 = -(2^(-4 + 3*\[CapitalDelta])*Gamma[1/2 + \[CapitalDelta]/2]*Gamma[\[CapitalDelta]/2]^3*Gamma[(1/2)*(-1 + \[CapitalDelta]) + \[CapitalDelta]\[Phi]]^2*Gamma[-(1/2) + 2*\[CapitalDelta]\[Phi]]*Gamma[-(1/2) + \[CapitalDelta] + 2*\[CapitalDelta]\[Phi]]*HypergeometricPFQ[{\[CapitalDelta]/2, \[CapitalDelta]/2, -(1/2) + \[CapitalDelta]/2 + \[CapitalDelta]\[Phi], -(1/2) + \[CapitalDelta]/2 + \[CapitalDelta]\[Phi], 1/4 + \[CapitalDelta]/2 + \[CapitalDelta]\[Phi], -1 + 2*\[CapitalDelta]\[Phi], -(3/2) + \[CapitalDelta] + 2*\[CapitalDelta]\[Phi]}, 
         {1/2 + \[CapitalDelta], -(3/4) + \[CapitalDelta]/2 + \[CapitalDelta]\[Phi], \[CapitalDelta]/2 + \[CapitalDelta]\[Phi], \[CapitalDelta]/2 + \[CapitalDelta]\[Phi], -(1/2) + \[CapitalDelta]/2 + 2*\[CapitalDelta]\[Phi], -(1/2) + \[CapitalDelta]/2 + 2*\[CapitalDelta]\[Phi]}, 1])/(Sqrt[Pi]*Gamma[2*\[CapitalDelta]]*Gamma[\[CapitalDelta]/2 + \[CapitalDelta]\[Phi]]^2*Gamma[-(1/2) + \[CapitalDelta]/2 + 2*\[CapitalDelta]\[Phi]]^2);
 b0[-1] = 0; b0[0] = \[Beta]0s0; 
    b0[(n_)?OddQ] := b0[n] = (R[(n - 1)/2] - (\[Rho][n - 2]*b0[n - 2] + \[Nu][n - 1]*b0[n - 1]))/\[Mu][n]; 
b0[(n_)?EvenQ] := b0[n] = -(\[Rho][n - 2]*b0[n - 2] + \[Nu][n - 1]*b0[n - 1])/\[Mu][n]; 
-\[Beta]SS[order]];


betafunctionalD2[ext_, order_, x_, OptionsPattern[]]:= Block[{
\[CapitalDelta]\[Phi] = SetPrecision[ext, OptionValue[PrecisionGoal]],
 \[CapitalDelta] = SetPrecision[x, OptionValue[PrecisionGoal]], \[Beta]0s0, b0, brule}, 
\[Beta]0s0 = -(2^(-4 + 3*\[CapitalDelta])*Gamma[1/2 + \[CapitalDelta]/2]*Gamma[\[CapitalDelta]/2]^3*Gamma[(1/2)*(-1 + \[CapitalDelta]) + \[CapitalDelta]\[Phi]]^2*Gamma[-(1/2) + 2*\[CapitalDelta]\[Phi]]*Gamma[-(1/2) + \[CapitalDelta] + 2*\[CapitalDelta]\[Phi]]*HypergeometricPFQ[{\[CapitalDelta]/2, \[CapitalDelta]/2, -(1/2) + \[CapitalDelta]/2 + \[CapitalDelta]\[Phi], -(1/2) + \[CapitalDelta]/2 + \[CapitalDelta]\[Phi], 1/4 + \[CapitalDelta]/2 + \[CapitalDelta]\[Phi], -1 + 2*\[CapitalDelta]\[Phi], -(3/2) + \[CapitalDelta] + 2*\[CapitalDelta]\[Phi]}, 
         {1/2 + \[CapitalDelta], -(3/4) + \[CapitalDelta]/2 + \[CapitalDelta]\[Phi], \[CapitalDelta]/2 + \[CapitalDelta]\[Phi], \[CapitalDelta]/2 + \[CapitalDelta]\[Phi], -(1/2) + \[CapitalDelta]/2 + 2*\[CapitalDelta]\[Phi], -(1/2) + \[CapitalDelta]/2 + 2*\[CapitalDelta]\[Phi]}, 1])/(Sqrt[Pi]*Gamma[2*\[CapitalDelta]]*Gamma[\[CapitalDelta]/2 + \[CapitalDelta]\[Phi]]^2*Gamma[-(1/2) + \[CapitalDelta]/2 + 2*\[CapitalDelta]\[Phi]]^2);
 b0[-1] = 0; b0[0] = \[Beta]0s0; 
    b0[(n_)?OddQ] := b0[n] = (R[(n - 1)/2] - (\[Rho][n - 2]*b0[n - 2] + \[Nu][n - 1]*b0[n - 1]))/\[Mu][n]; 
b0[(n_)?EvenQ] := b0[n] = -(\[Rho][n - 2]*b0[n - 2] + \[Nu][n - 1]*b0[n - 1])/\[Mu][n]; 
-\[Beta]SSD2[order]];



betafunctionallist[ext_, order_, x_,OptionsPattern[]]:= Block[{\[CapitalDelta]\[Phi] = SetPrecision[ext, OptionValue[PrecisionGoal]],
 \[CapitalDelta] = SetPrecision[x, OptionValue[PrecisionGoal]], \[Beta]0s0, b0, brule}, 
   \[Beta]0s0 = -(2^(-4 + 3*\[CapitalDelta])*Gamma[1/2 + \[CapitalDelta]/2]*Gamma[\[CapitalDelta]/2]^3*Gamma[(1/2)*(-1 + \[CapitalDelta]) + \[CapitalDelta]\[Phi]]^2*Gamma[-(1/2) + 2*\[CapitalDelta]\[Phi]]*Gamma[-(1/2) + \[CapitalDelta] + 2*\[CapitalDelta]\[Phi]]*HypergeometricPFQ[{\[CapitalDelta]/2, \[CapitalDelta]/2, -(1/2) + \[CapitalDelta]/2 + \[CapitalDelta]\[Phi], -(1/2) + \[CapitalDelta]/2 + \[CapitalDelta]\[Phi], 1/4 + \[CapitalDelta]/2 + \[CapitalDelta]\[Phi], -1 + 2*\[CapitalDelta]\[Phi], -(3/2) + \[CapitalDelta] + 2*\[CapitalDelta]\[Phi]}, 
         {1/2 + \[CapitalDelta], -(3/4) + \[CapitalDelta]/2 + \[CapitalDelta]\[Phi], \[CapitalDelta]/2 + \[CapitalDelta]\[Phi], \[CapitalDelta]/2 + \[CapitalDelta]\[Phi], -(1/2) + \[CapitalDelta]/2 + 2*\[CapitalDelta]\[Phi], -(1/2) + \[CapitalDelta]/2 + 2*\[CapitalDelta]\[Phi]}, 1])/(Sqrt[Pi]*Gamma[2*\[CapitalDelta]]*Gamma[\[CapitalDelta]/2 + \[CapitalDelta]\[Phi]]^2*Gamma[-(1/2) + \[CapitalDelta]/2 + 2*\[CapitalDelta]\[Phi]]^2);
 b0[-1] = 0; b0[0] = \[Beta]0s0; 
    b0[(n_)?OddQ] := b0[n] = (R[(n - 1)/2] - (\[Rho][n - 2]*b0[n - 2] + \[Nu][n - 1]*b0[n - 1]))/\[Mu][n]; 
b0[(n_)?EvenQ] := b0[n] = -(\[Rho][n - 2]*b0[n - 2] + \[Nu][n - 1]*b0[n - 1])/\[Mu][n]; 
-\[Beta]SS /@ Range[order]];


alphafunctional[ext_, order_, x_, OptionsPattern[]] := Block[{\[CapitalDelta]\[Phi] = SetPrecision[ext, 3/2*OptionValue[PrecisionGoal]],
 \[CapitalDelta] = SetPrecision[x, 3/2*OptionValue[PrecisionGoal]], \[Beta]0s0, b0, brule, f, \[Alpha]0s0pre, \[Alpha]0s0,a0, prec=OptionValue[PrecisionGoal]}, 
   f[s_] = (Gamma[s]^4*Gamma[-1/2+2*s+\[CapitalDelta]]*Gamma[s+\[CapitalDelta]/2-\[CapitalDelta]\[Phi]]^2*HypergeometricPFQ[{-1+2*s,1/4+s+\[CapitalDelta]/2,-3/2+2*s+\[CapitalDelta],s+\[CapitalDelta]/2-\[CapitalDelta]\[Phi],s+\[CapitalDelta]/2-\[CapitalDelta]\[Phi],-1/2+\[CapitalDelta]/2+\[CapitalDelta]\[Phi],-1/2+\[CapitalDelta]/2+\[CapitalDelta]\[Phi]},{-3/4+s+\[CapitalDelta]/2,1/2+\[CapitalDelta],2*s+\[CapitalDelta]/2-\[CapitalDelta]\[Phi],2*s+\[CapitalDelta]/2-\[CapitalDelta]\[Phi],-1/2+s+\[CapitalDelta]/2+\[CapitalDelta]\[Phi],-1/2+s+\[CapitalDelta]/2+\[CapitalDelta]\[Phi]},1])/(Gamma[2*s+\[CapitalDelta]/2-\[CapitalDelta]\[Phi]]^2*Gamma[-1/2+s+\[CapitalDelta]/2+\[CapitalDelta]\[Phi]]^2);
\[Alpha]0s0pre=(f[\[CapitalDelta]\[Phi]  + 1/2*10^(-prec/2)] - f[\[CapitalDelta]\[Phi]  - 1/2*10^(-prec/2)])/10^(-prec/2);
\[CapitalDelta]\[Phi] = SetPrecision[ext, prec];
 \[CapitalDelta] = SetPrecision[x, prec];
\[Beta]0s0 = -(2^(-4 + 3*\[CapitalDelta])*Gamma[1/2 + \[CapitalDelta]/2]*Gamma[\[CapitalDelta]/2]^3*Gamma[(1/2)*(-1 + \[CapitalDelta]) + \[CapitalDelta]\[Phi]]^2*Gamma[-(1/2) + 2*\[CapitalDelta]\[Phi]]*Gamma[-(1/2) + \[CapitalDelta] + 2*\[CapitalDelta]\[Phi]]*HypergeometricPFQ[{\[CapitalDelta]/2, \[CapitalDelta]/2, -(1/2) + \[CapitalDelta]/2 + \[CapitalDelta]\[Phi], -(1/2) + \[CapitalDelta]/2 + \[CapitalDelta]\[Phi], 1/4 + \[CapitalDelta]/2 + \[CapitalDelta]\[Phi], -1 + 2*\[CapitalDelta]\[Phi], -(3/2) + \[CapitalDelta] + 2*\[CapitalDelta]\[Phi]}, 
         {1/2 + \[CapitalDelta], -(3/4) + \[CapitalDelta]/2 + \[CapitalDelta]\[Phi], \[CapitalDelta]/2 + \[CapitalDelta]\[Phi], \[CapitalDelta]/2 + \[CapitalDelta]\[Phi], -(1/2) + \[CapitalDelta]/2 + 2*\[CapitalDelta]\[Phi], -(1/2) + \[CapitalDelta]/2 + 2*\[CapitalDelta]\[Phi]}, 1])/(Sqrt[Pi]*Gamma[2*\[CapitalDelta]]*Gamma[\[CapitalDelta]/2 + \[CapitalDelta]\[Phi]]^2*Gamma[-(1/2) + \[CapitalDelta]/2 + 2*\[CapitalDelta]\[Phi]]^2);
\[Alpha]0s0=-((2^(-5+3*\[CapitalDelta])*Gamma[1/2+\[CapitalDelta]/2]*Gamma[\[CapitalDelta]/2]*Gamma[(-1+\[CapitalDelta])/2+\[CapitalDelta]\[Phi]]^2*Gamma[-1/2+2*\[CapitalDelta]\[Phi]])/(Sqrt[Pi]*Gamma[2*\[CapitalDelta]]*Gamma[\[CapitalDelta]\[Phi]]^4)*\[Alpha]0s0pre-EulerGamma \[Beta]0s0);
	 b0[-1] = 0; b0[0] = \[Beta]0s0; 
    b0[(n_)?OddQ] := b0[n] = (R[(n - 1)/2] - (\[Rho][n - 2]*b0[n - 2] + \[Nu][n - 1]*b0[n - 1]))/\[Mu][n]; 
b0[(n_)?EvenQ] := b0[n] = -(\[Rho][n - 2]*b0[n - 2] + \[Nu][n - 1]*b0[n - 1])/\[Mu][n]; 
 a0[-1] = 0; a0[0] =\[Alpha]0s0; 
    a0[(n_)?OddQ] := a0[n] = (S[(n - 1)/2] - (\[Rho][n-2] a0[n-2] +\[Nu][n-1] a0[n-1]+\[Rho]p[n-2] b0[n-2]+\[Nu]p[n-1] b0[n-1]+\[Mu]p[n] b0[n]))/\[Mu][n]; 
a0[(n_)?EvenQ] := a0[n] = -(\[Rho][n-2] a0[n-2] +\[Nu][n-1] a0[n-1]+\[Rho]p[n-2] b0[n-2]+\[Nu]p[n-1] b0[n-1]+\[Mu]p[n] b0[n])/\[Mu][n]; 
-\[Alpha]SS@order];


alphafunctionallist[ext_, order_, x_, OptionsPattern[]] := Block[{\[CapitalDelta]\[Phi] = SetPrecision[ext, 3/2*OptionValue[PrecisionGoal]],
 \[CapitalDelta] = SetPrecision[x, 3/2*OptionValue[PrecisionGoal]], \[Beta]0s0, b0, brule, f, \[Alpha]0s0pre, \[Alpha]0s0,a0, prec=OptionValue[PrecisionGoal]}, 
  f[s_] = (Gamma[s]^4*Gamma[-1/2+2*s+\[CapitalDelta]]*Gamma[s+\[CapitalDelta]/2-\[CapitalDelta]\[Phi]]^2*HypergeometricPFQ[{-1+2*s,1/4+s+\[CapitalDelta]/2,-3/2+2*s+\[CapitalDelta],s+\[CapitalDelta]/2-\[CapitalDelta]\[Phi],s+\[CapitalDelta]/2-\[CapitalDelta]\[Phi],-1/2+\[CapitalDelta]/2+\[CapitalDelta]\[Phi],-1/2+\[CapitalDelta]/2+\[CapitalDelta]\[Phi]},{-3/4+s+\[CapitalDelta]/2,1/2+\[CapitalDelta],2*s+\[CapitalDelta]/2-\[CapitalDelta]\[Phi],2*s+\[CapitalDelta]/2-\[CapitalDelta]\[Phi],-1/2+s+\[CapitalDelta]/2+\[CapitalDelta]\[Phi],-1/2+s+\[CapitalDelta]/2+\[CapitalDelta]\[Phi]},1])/(Gamma[2*s+\[CapitalDelta]/2-\[CapitalDelta]\[Phi]]^2*Gamma[-1/2+s+\[CapitalDelta]/2+\[CapitalDelta]\[Phi]]^2);
\[Alpha]0s0pre=(f[\[CapitalDelta]\[Phi]  + 1/2*10^(-prec/2)] - f[\[CapitalDelta]\[Phi]  - 1/2*10^(-prec/2)])/10^(-prec/2);
\[CapitalDelta]\[Phi] = SetPrecision[ext, prec];
 \[CapitalDelta] = SetPrecision[x, prec];
\[Beta]0s0 = -(2^(-4 + 3*\[CapitalDelta])*Gamma[1/2 + \[CapitalDelta]/2]*Gamma[\[CapitalDelta]/2]^3*Gamma[(1/2)*(-1 + \[CapitalDelta]) + \[CapitalDelta]\[Phi]]^2*Gamma[-(1/2) + 2*\[CapitalDelta]\[Phi]]*Gamma[-(1/2) + \[CapitalDelta] + 2*\[CapitalDelta]\[Phi]]*HypergeometricPFQ[{\[CapitalDelta]/2, \[CapitalDelta]/2, -(1/2) + \[CapitalDelta]/2 + \[CapitalDelta]\[Phi], -(1/2) + \[CapitalDelta]/2 + \[CapitalDelta]\[Phi], 1/4 + \[CapitalDelta]/2 + \[CapitalDelta]\[Phi], -1 + 2*\[CapitalDelta]\[Phi], -(3/2) + \[CapitalDelta] + 2*\[CapitalDelta]\[Phi]}, 
         {1/2 + \[CapitalDelta], -(3/4) + \[CapitalDelta]/2 + \[CapitalDelta]\[Phi], \[CapitalDelta]/2 + \[CapitalDelta]\[Phi], \[CapitalDelta]/2 + \[CapitalDelta]\[Phi], -(1/2) + \[CapitalDelta]/2 + 2*\[CapitalDelta]\[Phi], -(1/2) + \[CapitalDelta]/2 + 2*\[CapitalDelta]\[Phi]}, 1])/(Sqrt[Pi]*Gamma[2*\[CapitalDelta]]*Gamma[\[CapitalDelta]/2 + \[CapitalDelta]\[Phi]]^2*Gamma[-(1/2) + \[CapitalDelta]/2 + 2*\[CapitalDelta]\[Phi]]^2);
\[Alpha]0s0=-((2^(-5+3*\[CapitalDelta])*Gamma[1/2+\[CapitalDelta]/2]*Gamma[\[CapitalDelta]/2]*Gamma[(-1+\[CapitalDelta])/2+\[CapitalDelta]\[Phi]]^2*Gamma[-1/2+2*\[CapitalDelta]\[Phi]])/(Sqrt[Pi]*Gamma[2*\[CapitalDelta]]*Gamma[\[CapitalDelta]\[Phi]]^4)*\[Alpha]0s0pre-EulerGamma \[Beta]0s0);
	 b0[-1] = 0; b0[0] = \[Beta]0s0; 
    b0[(n_)?OddQ] := b0[n] = (R[(n - 1)/2] - (\[Rho][n - 2]*b0[n - 2] + \[Nu][n - 1]*b0[n - 1]))/\[Mu][n]; 
b0[(n_)?EvenQ] := b0[n] = -(\[Rho][n - 2]*b0[n - 2] + \[Nu][n - 1]*b0[n - 1])/\[Mu][n]; 
 a0[-1] = 0; a0[0] =\[Alpha]0s0; 
    a0[(n_)?OddQ] := a0[n] = (S[(n - 1)/2] - (\[Rho][n-2] a0[n-2] +\[Nu][n-1] a0[n-1]+\[Rho]p[n-2] b0[n-2]+\[Nu]p[n-1] b0[n-1]+\[Mu]p[n] b0[n]))/\[Mu][n]; 
a0[(n_)?EvenQ] := a0[n] = -(\[Rho][n-2] a0[n-2] +\[Nu][n-1] a0[n-1]+\[Rho]p[n-2] b0[n-2]+\[Nu]p[n-1] b0[n-1]+\[Mu]p[n] b0[n])/\[Mu][n]; 
-\[Alpha]SS/@Range[0,order]];


functionallist[ext_, order_, x_, OptionsPattern[]] := Block[{\[CapitalDelta]\[Phi] = SetPrecision[ext, 3/2*OptionValue[PrecisionGoal]],
 \[CapitalDelta] = SetPrecision[x, 3/2*OptionValue[PrecisionGoal]], \[Beta]0s0, b0, brule, f, \[Alpha]0s0pre, \[Alpha]0s0,a0, prec=OptionValue[PrecisionGoal]}, 
   f[s_] = (Gamma[s]^4*Gamma[-1/2+2*s+\[CapitalDelta]]*Gamma[s+\[CapitalDelta]/2-\[CapitalDelta]\[Phi]]^2*HypergeometricPFQ[{-1+2*s,1/4+s+\[CapitalDelta]/2,-3/2+2*s+\[CapitalDelta],s+\[CapitalDelta]/2-\[CapitalDelta]\[Phi],s+\[CapitalDelta]/2-\[CapitalDelta]\[Phi],-1/2+\[CapitalDelta]/2+\[CapitalDelta]\[Phi],-1/2+\[CapitalDelta]/2+\[CapitalDelta]\[Phi]},{-3/4+s+\[CapitalDelta]/2,1/2+\[CapitalDelta],2*s+\[CapitalDelta]/2-\[CapitalDelta]\[Phi],2*s+\[CapitalDelta]/2-\[CapitalDelta]\[Phi],-1/2+s+\[CapitalDelta]/2+\[CapitalDelta]\[Phi],-1/2+s+\[CapitalDelta]/2+\[CapitalDelta]\[Phi]},1])/(Gamma[2*s+\[CapitalDelta]/2-\[CapitalDelta]\[Phi]]^2*Gamma[-1/2+s+\[CapitalDelta]/2+\[CapitalDelta]\[Phi]]^2);
\[Alpha]0s0pre=(f[\[CapitalDelta]\[Phi]  + 1/2*10^(-prec/2)] - f[\[CapitalDelta]\[Phi]  - 1/2*10^(-prec/2)])/10^(-prec/2);
\[CapitalDelta]\[Phi] = SetPrecision[ext, prec];
 \[CapitalDelta] = SetPrecision[x, prec];
\[Beta]0s0 = -(2^(-4 + 3*\[CapitalDelta])*Gamma[1/2 + \[CapitalDelta]/2]*Gamma[\[CapitalDelta]/2]^3*Gamma[(1/2)*(-1 + \[CapitalDelta]) + \[CapitalDelta]\[Phi]]^2*Gamma[-(1/2) + 2*\[CapitalDelta]\[Phi]]*Gamma[-(1/2) + \[CapitalDelta] + 2*\[CapitalDelta]\[Phi]]*HypergeometricPFQ[{\[CapitalDelta]/2, \[CapitalDelta]/2, -(1/2) + \[CapitalDelta]/2 + \[CapitalDelta]\[Phi], -(1/2) + \[CapitalDelta]/2 + \[CapitalDelta]\[Phi], 1/4 + \[CapitalDelta]/2 + \[CapitalDelta]\[Phi], -1 + 2*\[CapitalDelta]\[Phi], -(3/2) + \[CapitalDelta] + 2*\[CapitalDelta]\[Phi]}, 
         {1/2 + \[CapitalDelta], -(3/4) + \[CapitalDelta]/2 + \[CapitalDelta]\[Phi], \[CapitalDelta]/2 + \[CapitalDelta]\[Phi], \[CapitalDelta]/2 + \[CapitalDelta]\[Phi], -(1/2) + \[CapitalDelta]/2 + 2*\[CapitalDelta]\[Phi], -(1/2) + \[CapitalDelta]/2 + 2*\[CapitalDelta]\[Phi]}, 1])/(Sqrt[Pi]*Gamma[2*\[CapitalDelta]]*Gamma[\[CapitalDelta]/2 + \[CapitalDelta]\[Phi]]^2*Gamma[-(1/2) + \[CapitalDelta]/2 + 2*\[CapitalDelta]\[Phi]]^2);
\[Alpha]0s0=-((2^(-5+3*\[CapitalDelta])*Gamma[1/2+\[CapitalDelta]/2]*Gamma[\[CapitalDelta]/2]*Gamma[(-1+\[CapitalDelta])/2+\[CapitalDelta]\[Phi]]^2*Gamma[-1/2+2*\[CapitalDelta]\[Phi]])/(Sqrt[Pi]*Gamma[2*\[CapitalDelta]]*Gamma[\[CapitalDelta]\[Phi]]^4)*\[Alpha]0s0pre-EulerGamma \[Beta]0s0);
	 b0[-1] = 0; b0[0] = \[Beta]0s0; 
    b0[(n_)?OddQ] := b0[n] = (R[(n - 1)/2] - (\[Rho][n - 2]*b0[n - 2] + \[Nu][n - 1]*b0[n - 1]))/\[Mu][n]; 
b0[(n_)?EvenQ] := b0[n] = -(\[Rho][n - 2]*b0[n - 2] + \[Nu][n - 1]*b0[n - 1])/\[Mu][n]; 
 a0[-1] = 0; a0[0] =\[Alpha]0s0; 
    a0[(n_)?OddQ] := a0[n] = (S[(n - 1)/2] - (\[Rho][n-2] a0[n-2] +\[Nu][n-1] a0[n-1]+\[Rho]p[n-2] b0[n-2]+\[Nu]p[n-1] b0[n-1]+\[Mu]p[n] b0[n]))/\[Mu][n]; 
a0[(n_)?EvenQ] := a0[n] = -(\[Rho][n-2] a0[n-2] +\[Nu][n-1] a0[n-1]+\[Rho]p[n-2] b0[n-2]+\[Nu]p[n-1] b0[n-1]+\[Mu]p[n] b0[n])/\[Mu][n]; 
{-\[Beta]SS /@ Range[order],-\[Alpha]SS/@({0}~Join~Range[order])}];


functionallistPara[ext_, order_, xlist_, opt___]:=ParallelTable[Print[x];functionallist[ext,order,x],{x,xlist}];


functional[ext_, order_, x_, OptionsPattern[]] := Block[{\[CapitalDelta]\[Phi] = SetPrecision[ext, 3/2*OptionValue[PrecisionGoal]],
 \[CapitalDelta] = SetPrecision[x, 3/2*OptionValue[PrecisionGoal]], \[Beta]0s0, b0, brule, f, \[Alpha]0s0pre, \[Alpha]0s0,a0, prec=OptionValue[PrecisionGoal]}, 
   f[s_] = (Gamma[s]^4*Gamma[-1/2+2*s+\[CapitalDelta]]*Gamma[s+\[CapitalDelta]/2-\[CapitalDelta]\[Phi]]^2*HypergeometricPFQ[{-1+2*s,1/4+s+\[CapitalDelta]/2,-3/2+2*s+\[CapitalDelta],s+\[CapitalDelta]/2-\[CapitalDelta]\[Phi],s+\[CapitalDelta]/2-\[CapitalDelta]\[Phi],-1/2+\[CapitalDelta]/2+\[CapitalDelta]\[Phi],-1/2+\[CapitalDelta]/2+\[CapitalDelta]\[Phi]},{-3/4+s+\[CapitalDelta]/2,1/2+\[CapitalDelta],2*s+\[CapitalDelta]/2-\[CapitalDelta]\[Phi],2*s+\[CapitalDelta]/2-\[CapitalDelta]\[Phi],-1/2+s+\[CapitalDelta]/2+\[CapitalDelta]\[Phi],-1/2+s+\[CapitalDelta]/2+\[CapitalDelta]\[Phi]},1])/(Gamma[2*s+\[CapitalDelta]/2-\[CapitalDelta]\[Phi]]^2*Gamma[-1/2+s+\[CapitalDelta]/2+\[CapitalDelta]\[Phi]]^2);
\[Alpha]0s0pre=(f[\[CapitalDelta]\[Phi]  + 1/2*10^(-prec/2)] - f[\[CapitalDelta]\[Phi]  - 1/2*10^(-prec/2)])/10^(-prec/2);
\[CapitalDelta]\[Phi] = SetPrecision[ext, prec];
 \[CapitalDelta] = SetPrecision[x, prec];
\[Beta]0s0 = -(2^(-4 + 3*\[CapitalDelta])*Gamma[1/2 + \[CapitalDelta]/2]*Gamma[\[CapitalDelta]/2]^3*Gamma[(1/2)*(-1 + \[CapitalDelta]) + \[CapitalDelta]\[Phi]]^2*Gamma[-(1/2) + 2*\[CapitalDelta]\[Phi]]*Gamma[-(1/2) + \[CapitalDelta] + 2*\[CapitalDelta]\[Phi]]*HypergeometricPFQ[{\[CapitalDelta]/2, \[CapitalDelta]/2, -(1/2) + \[CapitalDelta]/2 + \[CapitalDelta]\[Phi], -(1/2) + \[CapitalDelta]/2 + \[CapitalDelta]\[Phi], 1/4 + \[CapitalDelta]/2 + \[CapitalDelta]\[Phi], -1 + 2*\[CapitalDelta]\[Phi], -(3/2) + \[CapitalDelta] + 2*\[CapitalDelta]\[Phi]}, 
         {1/2 + \[CapitalDelta], -(3/4) + \[CapitalDelta]/2 + \[CapitalDelta]\[Phi], \[CapitalDelta]/2 + \[CapitalDelta]\[Phi], \[CapitalDelta]/2 + \[CapitalDelta]\[Phi], -(1/2) + \[CapitalDelta]/2 + 2*\[CapitalDelta]\[Phi], -(1/2) + \[CapitalDelta]/2 + 2*\[CapitalDelta]\[Phi]}, 1])/(Sqrt[Pi]*Gamma[2*\[CapitalDelta]]*Gamma[\[CapitalDelta]/2 + \[CapitalDelta]\[Phi]]^2*Gamma[-(1/2) + \[CapitalDelta]/2 + 2*\[CapitalDelta]\[Phi]]^2);
\[Alpha]0s0=-((2^(-5+3*\[CapitalDelta])*Gamma[1/2+\[CapitalDelta]/2]*Gamma[\[CapitalDelta]/2]*Gamma[(-1+\[CapitalDelta])/2+\[CapitalDelta]\[Phi]]^2*Gamma[-1/2+2*\[CapitalDelta]\[Phi]])/(Sqrt[Pi]*Gamma[2*\[CapitalDelta]]*Gamma[\[CapitalDelta]\[Phi]]^4)*\[Alpha]0s0pre-EulerGamma \[Beta]0s0);
	 b0[-1] = 0; b0[0] = \[Beta]0s0; 
    b0[(n_)?OddQ] := b0[n] = (R[(n - 1)/2] - (\[Rho][n - 2]*b0[n - 2] + \[Nu][n - 1]*b0[n - 1]))/\[Mu][n]; 
b0[(n_)?EvenQ] := b0[n] = -(\[Rho][n - 2]*b0[n - 2] + \[Nu][n - 1]*b0[n - 1])/\[Mu][n]; 
 a0[-1] = 0; a0[0] =\[Alpha]0s0; 
    a0[(n_)?OddQ] := a0[n] = (S[(n - 1)/2] - (\[Rho][n-2] a0[n-2] +\[Nu][n-1] a0[n-1]+\[Rho]p[n-2] b0[n-2]+\[Nu]p[n-1] b0[n-1]+\[Mu]p[n] b0[n]))/\[Mu][n]; 
a0[(n_)?EvenQ] := a0[n] = -(\[Rho][n-2] a0[n-2] +\[Nu][n-1] a0[n-1]+\[Rho]p[n-2] b0[n-2]+\[Nu]p[n-1] b0[n-1]+\[Mu]p[n] b0[n])/\[Mu][n]; 
{-\[Beta]SS[order],-\[Alpha]SS[order]}];

betafunctionalderi2[ext_, order_, m_, OptionsPattern[]]:= Block[{\[CapitalDelta]\[Phi] = SetPrecision[ext, OptionValue[PrecisionGoal]],
 \[CapitalDelta] = SetPrecision[2 ext+2 m, OptionValue[PrecisionGoal]], \[Beta]0s0, b0,\[Rho],\[Nu],\[Mu],R,\[Beta]SS,bs0,bCS,A},
 \[Beta]0s0 = -(2^(-4 + 3*\[CapitalDelta])*Gamma[1/2 + \[CapitalDelta]/2]*Gamma[\[CapitalDelta]/2]^3*Gamma[(1/2)*(-1 + \[CapitalDelta]) + \[CapitalDelta]\[Phi]]^2*Gamma[-(1/2) + 2*\[CapitalDelta]\[Phi]]*Gamma[-(1/2) + \[CapitalDelta] + 2*\[CapitalDelta]\[Phi]]*HypergeometricPFQ[{\[CapitalDelta]/2, \[CapitalDelta]/2, -(1/2) + \[CapitalDelta]/2 + \[CapitalDelta]\[Phi], -(1/2) + \[CapitalDelta]/2 + \[CapitalDelta]\[Phi], 1/4 + \[CapitalDelta]/2 + \[CapitalDelta]\[Phi], -1 + 2*\[CapitalDelta]\[Phi], -(3/2) + \[CapitalDelta] + 2*\[CapitalDelta]\[Phi]}, 
         {1/2 + \[CapitalDelta], -(3/4) + \[CapitalDelta]/2 + \[CapitalDelta]\[Phi], \[CapitalDelta]/2 + \[CapitalDelta]\[Phi], \[CapitalDelta]/2 + \[CapitalDelta]\[Phi], -(1/2) + \[CapitalDelta]/2 + 2*\[CapitalDelta]\[Phi], -(1/2) + \[CapitalDelta]/2 + 2*\[CapitalDelta]\[Phi]}, 1])/(Sqrt[Pi]*Gamma[2*\[CapitalDelta]]*Gamma[\[CapitalDelta]/2 + \[CapitalDelta]\[Phi]]^2*Gamma[-(1/2) + \[CapitalDelta]/2 + 2*\[CapitalDelta]\[Phi]]^2);
 b0[-1] = 0; b0[0] = \[Beta]0s0; 
R[n_]:=(\[Pi]^(1/2) (-1)^(1-2n) Gamma[n+\[CapitalDelta]\[Phi]]^4 Gamma[-1/2+n+2\[CapitalDelta]\[Phi]]^2)/((n!)^2 Gamma[\[CapitalDelta]\[Phi]]^4 Gamma[2(n+\[CapitalDelta]\[Phi])]Gamma[2(n+\[CapitalDelta]\[Phi])-1/2]);
\[Mu][n_]:=-n^2;
\[Nu][n_]:=(-1+\[CapitalDelta]) \[CapitalDelta]+\[CapitalDelta]\[Phi]+1/2 n (-1+n+4 \[CapitalDelta]\[Phi]);
\[Rho][n_]:=-(((n+2 \[CapitalDelta]\[Phi])^2 (-1+n+4 \[CapitalDelta]\[Phi])^2)/(4 (-1+2 n+4 \[CapitalDelta]\[Phi]) (1+2 n+4 \[CapitalDelta]\[Phi])));
bs0[n_]=((-1)^(2*n)*Sqrt[Pi]*Gamma[n+\[CapitalDelta]\[Phi]]^4*Gamma[-1/2+n+2*\[CapitalDelta]\[Phi]]^2)/((2*n-\[CapitalDelta]+2*\[CapitalDelta]\[Phi])*(-1+2*n+\[CapitalDelta]+2*\[CapitalDelta]\[Phi])*n!^2*Gamma[\[CapitalDelta]\[Phi]]^4*Gamma[2*(n+\[CapitalDelta]\[Phi])]*Gamma[-1/2+2*n+2*\[CapitalDelta]\[Phi]]);
bCS[\[CapitalDelta]\[Phi]_,n_]=(2^(-1-2 (n+\[CapitalDelta]+\[CapitalDelta]\[Phi])) Sqrt[\[Pi]] Gamma[2 \[CapitalDelta]] Gamma[n+\[CapitalDelta]\[Phi]]^3 ((Sqrt[\[Pi]] (\[CapitalDelta]-2 \[CapitalDelta]\[Phi]) (-1+\[CapitalDelta]+2 \[CapitalDelta]\[Phi]) Gamma[\[CapitalDelta]/2-\[CapitalDelta]\[Phi]]^2)/Gamma[(1+\[CapitalDelta])/2+\[CapitalDelta]\[Phi]]^2-(32 \[Beta]0s0 Gamma[1+\[CapitalDelta]/2-\[CapitalDelta]\[Phi]]^2 Gamma[2 \[CapitalDelta]\[Phi]])/(Gamma[1/2 (-1+\[CapitalDelta])+\[CapitalDelta]\[Phi]]^2 Gamma[-(1/2)+2 \[CapitalDelta]\[Phi]])) Gamma[-(1/2)+n+2 \[CapitalDelta]\[Phi]]^2)/(Gamma[1+n]^2 Gamma[\[CapitalDelta]/2]^4 Gamma[1/2+n+\[CapitalDelta]\[Phi]] Gamma[-(1/2)+2 n+2 \[CapitalDelta]\[Phi]]);
A[\[CapitalDelta]_,\[CapitalDelta]\[Phi]_]=2/(Pi^2)*(4^(-2+\[CapitalDelta])*Gamma[\[CapitalDelta]/2]^4*Gamma[(-1+\[CapitalDelta])/2+\[CapitalDelta]\[Phi]]^2*((\[Pi]^2 )/Gamma[1/2 (2+\[CapitalDelta]-2 \[CapitalDelta]\[Phi])]^2))/(Gamma[2*\[CapitalDelta]]*Gamma[\[CapitalDelta]\[Phi]]^4);
\[Beta]SS[n_]= bs0[n]/A[\[CapitalDelta],\[CapitalDelta]\[Phi]]+2* b0[2n]/A[\[CapitalDelta],\[CapitalDelta]\[Phi]]+2bCS[\[CapitalDelta]\[Phi],n];
 b0[(n_)?OddQ] := b0[n] = (R[(n - 1)/2] - (\[Rho][n - 2]*b0[n - 2] + \[Nu][n - 1]*b0[n - 1]))/\[Mu][n];
b0[(n_)?EvenQ] := b0[n] = -(\[Rho][n - 2]*b0[n - 2] + \[Nu][n - 1]*b0[n - 1])/\[Mu][n]; 
-\[Beta]SS[order]];

tilde0D2[a_,\[CapitalDelta]_]:=2^(-5-2 \[CapitalDelta]) Sqrt[\[Pi]]Gamma[-1+a]^3 Gamma[-(5/2)+2 a] (-(1/(Gamma[1/2 (-2+\[CapitalDelta])]^4)) 2048 Gamma[1-a+\[CapitalDelta]/2]^2 Gamma[2 (-2+\[CapitalDelta])] ((4^(2-a) Sqrt[\[Pi]] (2+1/(-10 a+4 a^2-(-5+\[CapitalDelta]) \[CapitalDelta])))/(Gamma[-(1/2)+a] Gamma[-(5/2)+a+\[CapitalDelta]/2]^2)-1/2 (9+2 a-10 \[CapitalDelta]+2 \[CapitalDelta]^2) Gamma[-1+a] Gamma[-(11/4)+a+\[CapitalDelta]/2] Gamma[1/2 (-2+\[CapitalDelta])]^2 Gamma[-(9/2)+2 a+\[CapitalDelta]] HypergeometricPFQRegularized[{1/2 (-2+\[CapitalDelta]),1/2 (-2+\[CapitalDelta]),-3+2 a,1/2 (-5+2 a+\[CapitalDelta]),1/2 (-5+2 a+\[CapitalDelta]),-(7/4)+a+\[CapitalDelta]/2,-(11/2)+2 a+\[CapitalDelta]},{-(3/2)+\[CapitalDelta],-(11/4)+a+\[CapitalDelta]/2,-2+a+\[CapitalDelta]/2,-2+a+\[CapitalDelta]/2,1/2 (-7+4 a+\[CapitalDelta]),1/2 (-7+4 a+\[CapitalDelta])},1])+(1/((-3+2 \[CapitalDelta]) (1+2 \[CapitalDelta]) Gamma[\[CapitalDelta]/2]^4)) 64 (-1-\[CapitalDelta]+\[CapitalDelta]^2) Gamma[2-a+\[CapitalDelta]/2]^2 Gamma[2 \[CapitalDelta]] ((4^(2-a) Sqrt[\[Pi]] (2+1/(6-10 a+4 a^2+\[CapitalDelta]-\[CapitalDelta]^2)))/(Gamma[-(1/2)+a] Gamma[-(3/2)+a+\[CapitalDelta]/2]^2)-1/2 (-3+2 a-2 \[CapitalDelta]+2 \[CapitalDelta]^2) Gamma[-1+a] Gamma[-(7/4)+a+\[CapitalDelta]/2] Gamma[\[CapitalDelta]/2]^2 Gamma[-(5/2)+2 a+\[CapitalDelta]] HypergeometricPFQRegularized[{\[CapitalDelta]/2,\[CapitalDelta]/2,-3+2 a,1/2 (-3+2 a+\[CapitalDelta]),1/2 (-3+2 a+\[CapitalDelta]),-(3/4)+a+\[CapitalDelta]/2,-(7/2)+2 a+\[CapitalDelta]},{1/2+\[CapitalDelta],-(7/4)+a+\[CapitalDelta]/2,-1+a+\[CapitalDelta]/2,-1+a+\[CapitalDelta]/2,1/2 (-5+4 a+\[CapitalDelta]),1/2 (-5+4 a+\[CapitalDelta])},1])-(\[CapitalDelta]^2 (1+\[CapitalDelta])^2 Gamma[3-a+\[CapitalDelta]/2]^2 Gamma[2 (2+\[CapitalDelta])] ((4^(2-a) Sqrt[\[Pi]] (2+1/(4-10 a+4 a^2-3 \[CapitalDelta]-\[CapitalDelta]^2)))/(Gamma[-(1/2)+a] Gamma[-(1/2)+a+\[CapitalDelta]/2]^2)-1/2 (1+2 a+6 \[CapitalDelta]+2 \[CapitalDelta]^2) Gamma[-1+a] Gamma[-(3/4)+a+\[CapitalDelta]/2] Gamma[(2+\[CapitalDelta])/2]^2 Gamma[-(1/2)+2 a+\[CapitalDelta]] HypergeometricPFQRegularized[{(2+\[CapitalDelta])/2,(2+\[CapitalDelta])/2,-3+2 a,1/2 (-1+2 a+\[CapitalDelta]),1/2 (-1+2 a+\[CapitalDelta]),1/4+a+\[CapitalDelta]/2,-(3/2)+2 a+\[CapitalDelta]},{5/2+\[CapitalDelta],-(3/4)+a+\[CapitalDelta]/2,a+\[CapitalDelta]/2,a+\[CapitalDelta]/2,1/2 (-3+4 a+\[CapitalDelta]),1/2 (-3+4 a+\[CapitalDelta])},1]))/(2 (1+2 \[CapitalDelta])^2 (-3+4 \[CapitalDelta]+4 \[CapitalDelta]^2) Gamma[(2+\[CapitalDelta])/2]^4));
tilde0D2N[a_,n_Integer,OptionsPattern[]]:=tilde0D2[a,2 a+SetPrecision[2,OptionValue[PrecisionGoal]]*n];

cb1[a_][\[CapitalDelta]_,z_]:=z^(\[CapitalDelta]-2a) Hypergeometric2F1[\[CapitalDelta],\[CapitalDelta],2\[CapitalDelta],z];
afermion[a_][n_]:=(2Pochhammer[2a,2n+1]^2)/((2n+1)!Pochhammer[4a+2n,2n+1]);
am[a_][m_]:=(2 Pochhammer[2 a,2 m]^2)/((2 m)! Pochhammer[-1+4 a+2 m,2 m]);
gfermion[a_][z_]:=1/z^(2a)+1/(1-z)^(2a)-1;
afree[a_][\[CapitalDelta]_]:=(2 Gamma[\[CapitalDelta]]^2 Gamma[\[CapitalDelta]+2a-1])/(Gamma[2a]^2 Gamma[2\[CapitalDelta]-1]Gamma[\[CapitalDelta]-2a+1]);
dn[a_][n_]:=(Pochhammer[a,n]^4*Pochhammer[4a-1,2n])/(n!^2Pochhammer[2a,n]^2Pochhammer[4a+2n-1,2n]);
cn[a_][n_]:=-((2^(-3+6*a+4*n)*Gamma[1/2+a]*Gamma[a+n]^4*Gamma[-1/2+2*a+n]^2*(HarmonicNumber[n]-2*HarmonicNumber[-1+a+n]+HarmonicNumber[-1+2*a+n]-2*HarmonicNumber[-2+4*a+2*n]+2*HarmonicNumber[-2+4*a+4*n]))/(Pi*n!^2*Gamma[a]^3*Gamma[-1/2+2*a]*Gamma[-1+4*a+4*n]));
betatildefunctional[ext_, order_, x_, opts___]:=betafunctional[ext, order, x, opts]+dn[ext][order]betatildezero[ext,x, opts];
alphatildefunctional[ext_, order_, x_, opts___]:=alphafunctional[ext, order, x, opts]+cn[ext][order]betatildezero[ext,x, opts];

End[]
EndPackage[]

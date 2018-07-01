# RT-model-for-a-vertical-inhomogeneous-medium
the code for the paper named "Bidirectional radiative intensity calculation based on Eddington approximation in a vertically inhomogeneous medium"

# In the code

Input:

"n" is number of the level, "n=2" mean the code is suit for 1 lay (2 levels).

"u0" is cosine of solar zentih angel.

"u" is cosine of zentih angel.

"Del_varph" is difference between the azimuth angel and the the solar azimuth angel

"t0" is optical depth of the medium.

"ww" is single scattering albedo at "t0/2." in the medium.

"g" is asymmetry factor at "t0/2" in the medium.

"sigmaw"/"sigmag" are small coefficients in Eq. (4).

Output:

"fuu(n)" is upward flux of the vertically inhomogeneous medium at level n. (n=1,2)

"fdu(n)" is downward flux of the vertically inhomogeneous medium at level n. (n=1,2)

"TI" is azimuthally averaged radiative intensity at top of the inhomogeneous medium.

"TI_AZI" is azimuthally independent radiative intensity at top of the inhomogeneous medium.


module AdaptiveRejectionSampling

open System

open System

#if INTERACTIVE
open RDotNet
open RProvider
open RProvider.graphics
open RProvider.grDevices
#endif

/// Squeezing function, forms lower bound on the function
type LowerHull = 
    { M : float
      B : float
      Left : float
      Right: float }

/// Envelope function, forms upper bound on the function
type UpperHull =
    { M : float
      B : float
      Left : float
      Right: float
      Pr : float }

/// Compute hulls from a set of points and their function values
let arsComputeHulls domain (S: float[]) (fS: float[]) =
    // lower piecewise-linear hull
    let lowerHull = 
        Seq.zip S fS
        |> Seq.pairwise
        |> Seq.map (fun ((s1, fs1), (s2, fs2)) ->
            let m = (fs2 - fs1)/(s2 - s1)
            let b = fs1 - m * s1
            {M = m; B = b; Left = s1; Right = s2}:LowerHull)
        |> Array.ofSeq 

    // upper hull:
    // first line from the domain boundary
    let upperHullBoundaryBeg =
        let m = (fS.[1] - fS.[0])/(S.[1] - S.[0])
        let b = fS.[0] - m*S.[0]
        let pr = exp(b)/m * (exp(m*S.[0]) - 0.0)  // integrating from -infinity (or boundary point?)
        [| { M = m; B = b; Pr = pr; Left = fst domain; Right = S.[0] } |]

    // upper hull
    let upperHullBeg =
        let m = (fS.[2] - fS.[1])/(S.[2] - S.[1])
        let b = fS.[1] - m * S.[1]
        let pr = exp(b)/m * (exp(m * S.[1]) - exp(m * S.[0]))
        [| { M = m; B = b; Pr = pr; Left = S.[0]; Right = S.[1] } |]

    // interior lines
    // there are two lines between each abscissa
    let upperHullMid = 
      [|  for li in 1..S.Length-3 do
            let m1 = (fS.[li] - fS.[li-1])/(S.[li] - S.[li-1])
            let b1 = fS.[li] - m1 * S.[li]

            let m2 = (fS.[li+2]-fS.[li+1])/(S.[li+2] - S.[li+1])
            let b2 = fS.[li+1] - m2*S.[li+1]

            let ix = (b1 - b2)/(m2 - m1)   // intersection of the two lines

            let pr1 = exp(b1)/m1 * (exp(m1 * ix) - exp(m1*S.[li]))
            yield { M = m1; B = b1; Pr = pr1; Left = S.[li]; Right = ix }

            let pr2 = exp(b2)/m2 * ( exp(m2 * S.[li+1]) - exp(m2 * ix) )
            yield {M = m2; B = b2; Pr = pr2; Left = ix; Right = S.[li + 1] }
            |]

    // second last line
    let upperHullEnd =
        let m = (fS.[fS.Length-2] - fS.[fS.Length - 3])/(S.[fS.Length - 2] - S.[fS.Length - 3])
        let b = fS.[fS.Length - 2] - m*S.[fS.Length - 2]
        let pr = exp(b)/m * ( exp(m*S.[fS.Length-1]) - exp(m*S.[fS.Length-2]))
        [| { M = m; B = b; Pr = pr; Left = S.[fS.Length - 2]; Right = S.[fS.Length - 1] } |]

    // last line to the end of the domain
    let upperHullBoundaryEnd =
        let m = (fS.[fS.Length - 1] - fS.[fS.Length - 2])/(S.[fS.Length - 1] - S.[fS.Length - 2])
        let b = fS.[fS.Length - 1] - m * S.[fS.Length - 1]
        let pr = exp(b)/m * (0.0 - exp(m * S.[fS.Length - 1]))
        [| {M = m; B = b; Pr = pr; Left = S.[fS.Length - 1]; Right = snd domain } |]

    let upperHullUnnorm = 
        [| upperHullBoundaryBeg; 
           upperHullBeg; 
           upperHullMid; 
           upperHullEnd; 
           upperHullBoundaryEnd |]
         |> Array.concat 

    let Z = upperHullUnnorm |> Array.sumBy (fun uh -> uh.Pr)
    
    let upperHull = 
        upperHullUnnorm 
        |> Array.map (fun uh -> { uh with Pr = uh.Pr/Z })

    lowerHull, upperHull

let arsSampleUpperHull (upperHull: UpperHull[]) (rnd:Random) =
    let cdf = 
        Array.init upperHull.Length (fun i -> 
            upperHull.[..i] |> Array.map (fun uh -> uh.Pr) |> Array.sum)

    // randomly choose a line segment
    let U = rnd.NextDouble()
    let li = cdf |> Array.findIndex (fun s -> U < s)    // index of a line segment

    // sample along that line segment
    let U' = rnd.NextDouble()

    let m = upperHull.[li].M
    let b = upperHull.[li].B
    let left = upperHull.[li].Left
    let right = upperHull.[li].Right

    let x = log (U' * (exp(m*right) - exp(m*left)) + exp(m*left)) / m

    if Double.IsInfinity(x) || Double.IsNaN(x) then
        failwith "Sampled an infinite or NaN x."
    x


let (|OutsideOnLeft|OutsideOnRight|InsideInterval|) (lowerHull:LowerHull[], x) =
    if x < (lowerHull |> Array.map (fun lh -> lh.Left) |> Array.min)
    then OutsideOnLeft () 
    elif x > (lowerHull |> Array.map (fun lh -> lh.Right) |> Array.max)
    then OutsideOnRight () 
    else
        let intervalIdx =      
            lowerHull 
            |> Array.findIndex (fun h -> h.Left <= x && x <= h.Right)    
        InsideInterval (intervalIdx)  

/// Evaluate hull values at x
let arsEvalHulls (x:float) (lowerHull:LowerHull[]) (upperHull:UpperHull[]) =
    // lower bound
    let lhVal =
        match (lowerHull, x) with
        | OutsideOnLeft -> -infinity
        | OutsideOnRight -> -infinity
        | InsideInterval intervalIdx ->
            let lh = lowerHull.[intervalIdx]
            lh.M * x + lh.B
    
    // upper bound
    let uhVal =
        let uh = upperHull |> Array.find (fun h -> x >= h.Left && x <= h.Right)
        uh.M * x + uh.B

    lhVal, uhVal

/// Plot the upper and lower hulls, for debugging purposes only.
/// Requires RProvider and RDotNet libraries.
let arsPlot (upperHull:UpperHull[]) (lowerHull:LowerHull[]) domain (S:float[]) (fS:float[]) func =
    let Swidth = S.[S.Length - 1] - S.[0]
    let plotStep = Swidth/1000.0
    // plot this much before a and past b, if the domain is infinite
    let ext = 0.15 * Swidth

    let left = 
        if Double.IsNegativeInfinity(fst domain) 
        then S.[0] - ext 
        else S.[0]

    let right = 
        if Double.IsPositiveInfinity(snd domain) 
        then S.[S.Length-1] + ext
        else S.[S.Length-1]

    let x = [| left .. plotStep .. right |]
    let fx = x |> Array.map func

    #if INTERACTIVE
    R.x11() |> ignore
    namedParams ["x", box x; "y", box fx; "type", box "l"] |> R.plot |> ignore
    namedParams ["x", box S; "y", box fS] |> R.points   |> ignore
    #endif

    // Plot lower hull
    for li in 0..S.Length-2 do
        let m = lowerHull.[li].M
        let b = lowerHull.[li].B
        let x = [| lowerHull.[li].Left .. plotStep .. lowerHull.[li].Right |]
        let y = x |> Array.map (fun xi -> xi * m + b)
        #if INTERACTIVE
        namedParams ["x", box x; "y", box y; "col", box "blue"] 
        |> R.lines |> ignore
        #endif
        ()

    // Plot upper hull
    if Double.IsPositiveInfinity(fst domain) then
        let x = [| upperHull.[0].Right - ext .. plotStep .. upperHull.[0].Right |]
        let m = upperHull.[0].M
        let b = upperHull.[0].B
        let y = x |> Array.map (fun xi -> xi * m + b)
        #if INTERACTIVE
        namedParams ["x", box x; "y", box y; "col", box "red"]
        |> R.lines |> ignore
        #endif 
        ()

    // middle lines
    for li in 1..upperHull.Length-2 do
        let x = [| upperHull.[li].Left .. plotStep .. upperHull.[li].Right |]
        let m = upperHull.[li].M
        let b = upperHull.[li].B
        let y = x |> Array.map (fun xi -> xi * m + b)
        #if INTERACTIVE
        namedParams ["x", box x; "y", box y; "col", box "red"]
        |> R.lines |> ignore
        #endif
        ()

    // last line
    if Double.IsPositiveInfinity(snd domain) then
        let x = [| upperHull.[upperHull.Length - 1].Left .. plotStep .. upperHull.[upperHull.Length - 1].Left + ext|]
        let m = upperHull.[upperHull.Length - 1].M
        let b = upperHull.[upperHull.Length - 1].B
        let y = x |> Array.map (fun xi -> xi * m + b)
        #if INTERACTIVE
        namedParams ["x", box x; "y", box y; "col", box "red"]
        |> R.lines |> ignore
        #endif
        ()


// ================================================================================================
/// Adaptive rejection sampling (derivative-free). Based on matlab implementation from 
/// [PMTK library](https://github.com/probml).
/// 
/// # Parameters
///
/// * `func` - Function which computes logarithm of a likelihood function
/// * `a`, `b`  - Starting points, a < b. 
/// * `domain` - Tuple of `func` boundaries. If the left boundary is -infinity, 
///    derivative of `f` at `a` must be positive. If the right boundary is infinity,
///    derivative of `f` at `b` must be negative. 
/// * `nSamples` - Total number of samples that the sampler should return.
///
let adaptiveRejectionSampling func a b domain nSamples =

    let rnd = System.Random()

    let numDerivStep = 1e-3

    // check concavity and locations of initial points
    let a', b' = a + numDerivStep, b - numDerivStep
    let begGradientPositive = func a' - func a >= 0.0
    let endGradientNegative = func b - func b' <= 0.0

    if b <= a || (fst domain = -infinity && not begGradientPositive)
              || (snd domain = infinity && not endGradientNegative)
    then failwith "Incorrect initial points"

    // initialize mesh on which to create upper and lower hulls
    let nInitialMeshPoints = 3
    let S = 
        [| [|a|]; [|a' .. (b'-a')/(float nInitialMeshPoints + 1.0) .. b' |]; [|b|]|]
        |> Array.concat
    let fS = S |> Array.map func

    let lowerHull, upperHull = arsComputeHulls domain S fS

    //arsPlot upperHull lowerHull domain S fS func

    let rec getSamples lowerHull upperHull S (samples: float list) =
        if samples.Length = nSamples then 
            samples
        else
            // sample x from Hull
            let x = arsSampleUpperHull upperHull rnd
            let lhVal, uhVal = arsEvalHulls x lowerHull upperHull 

            let U = rnd.NextDouble()

            // three cases for acception/rejection
            let newSamples, meshChanged =
                if U <= exp(lhVal - uhVal) then
                    // accept, U is below lower bound
                    x :: samples, false
                elif U <= exp(func(x) - uhVal) then
                    // accept, u is between lower bound and f
                    x :: samples, true
                else
                    // reject, u is between f and upper bound
                    samples, true

            // recompute hulls if changed
            let newS, newLowerHull, newUpperHull =
                if meshChanged then
                    // recompute hulls
                    let nS = Array.append [|x|] S |> Array.sort
                    let nfS = nS |> Array.map func
                    let lHs, uHs = arsComputeHulls domain nS nfS
                    nS, lHs, uHs
                else
                    S, lowerHull, upperHull
                       
            getSamples newLowerHull newUpperHull newS newSamples 
                
    getSamples lowerHull upperHull S [] 


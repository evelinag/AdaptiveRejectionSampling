module AdaptiveRejectionSampling.ArsRewrite

open System
#if INTERACTIVE
#I "packages/RProvider.1.1.6"
#load "RProvider.fsx"
#endif

open System
open RDotNet
open RProvider
open RProvider.graphics
open RProvider.grDevices


type LowerHull = 
    { M : float
      B : float
      Left : float
      Right: float }

type UpperHull =
    { M : float
      B : float
      Left : float
      Right: float
      Pr : float }

let arsComputeHulls domain (S: float[]) (fS: float[]) =

    // lower piecewise-linear hull
    let (lowerHull: LowerHull[]) = 
        [| for li in 0..S.Length - 2 ->
            let m = (fS.[li + 1] - fS.[li])/(S.[li + 1] - S.[li])
            let b = fS.[li] - m * S.[li]
            { M = m
              B = b
              Left = S.[li]
              Right = S.[li+1]
            } |]

    // first line from the domain boundary
    let upperHullBoundaryBeg =
        let m = (fS.[1] - fS.[0])/(S.[1] - S.[0])
        let b = fS.[0] - m*S.[0]
        let pr = exp(b)/m * (exp(m*S.[0]) - 0.0)  // integrating from -infinity (or boundary point?)
        [| { M = m; B = b; Pr = pr; Left = fst domain; Right = S.[0] } |]

    // upper hull
    let m0 = (fS.[2] - fS.[1])/(S.[2] - S.[1])
    let b0 = fS.[1] - m0 * S.[1]
    let pr0 = exp(b0)/m0 * (exp(m0*S.[1]) - exp(m0*S.[0]))
    let upperHullBeg = 
        [| {M = m0; B = b0; Pr = pr0; Left = S.[0]; Right = S.[1] } |]

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

            let i2 = (li - 1)*2 + 1
            let pr2 = exp(b2)/m2 * ( exp(m2 * S.[li+1]) - exp(m2 * ix) )
            yield {M = m2; B = b2; Pr = pr2; Left = ix; Right = S.[li + 1] }
            |]

    // second last line
    let m = (fS.[fS.Length-2] - fS.[fS.Length - 3])/(S.[fS.Length - 2] - S.[fS.Length - 3])
    let b = fS.[fS.Length - 2] - m*S.[fS.Length - 2]
    let pr = exp(b)/m * ( exp(m*S.[fS.Length-1]) - exp(m*S.[fS.Length-2]))
    let upperHullEnd = 
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

/// Evaluate hull values at x
let arsEvalHulls (x:float) (lowerHull:LowerHull[]) (upperHull:UpperHull[]) =
    // lower bound
    let lhVal =
        if x < (lowerHull |> Array.map (fun lh -> lh.Left) |> Array.min) then
            Double.NegativeInfinity
        elif x > (lowerHull |> Array.map (fun lh -> lh.Right) |> Array.max) then
            Double.NegativeInfinity
        else
            let lh = 
                lowerHull 
                |> Array.find (fun h -> h.Left <= x && x <= h.Right)
            lh.M * x + lh.B
    
    // upper bound
    let uhVal =
        let uh = upperHull |> Array.find (fun h -> x >= h.Left && x <= h.Right)
        uh.M * x + uh.B

    lhVal, uhVal


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

    R.x11() |> ignore
    namedParams ["x", box x; "y", box fx; "type", box "l"] |> R.plot |> ignore
    namedParams ["x", box S; "y", box fS] |> R.points   |> ignore

    // Plot lower hull
    for li in 0..S.Length-2 do
        let m = lowerHull.[li].M
        let b = lowerHull.[li].B
        let x = [| lowerHull.[li].Left .. plotStep .. lowerHull.[li].Right |]
        let y = x |> Array.map (fun xi -> xi * m + b)
        namedParams ["x", box x; "y", box y; "col", box "blue"] 
        |> R.lines |> ignore

    // Plot upper hull
    if Double.IsPositiveInfinity(fst domain) then
        let x = [| upperHull.[0].Right - ext .. plotStep .. upperHull.[0].Right |]
        let m = upperHull.[0].M
        let b = upperHull.[0].B
        let y = x |> Array.map (fun xi -> xi * m + b)
        namedParams ["x", box x; "y", box y; "col", box "red"]
        |> R.lines |> ignore

    // middle lines
    for li in 1..upperHull.Length-2 do
        let x = [| upperHull.[li].Left .. plotStep .. upperHull.[li].Right |]
        let m = upperHull.[li].M
        let b = upperHull.[li].B
        let y = x |> Array.map (fun xi -> xi * m + b)
        namedParams ["x", box x; "y", box y; "col", box "red"]
        |> R.lines |> ignore

    // last line
    if Double.IsPositiveInfinity(snd domain) then
        let x = [| upperHull.[upperHull.Length - 1].Left .. plotStep .. upperHull.[upperHull.Length - 1].Left + ext|]
        let m = upperHull.[upperHull.Length - 1].M
        let b = upperHull.[upperHull.Length - 1].B
        let y = x |> Array.map (fun xi -> xi * m + b)
        namedParams ["x", box x; "y", box y; "col", box "red"]
        |> R.lines |> ignore

// ================================================================================================
// ================================================================================================
let ars func a b domain nSamples =

    let rnd = System.Random()

    let numDerivStep = 1e-3

    // initialize mesh on which to create upper and lower hulls
    let nInitialMeshPoints = 3
    let S = 
        let a', b' = a + numDerivStep, b - numDerivStep
        [| [|a|]; [|a' .. (b'-a')/(float nInitialMeshPoints + 1.0) .. b' |]; [|b|]|]
        |> Array.concat
    let fS = S |> Array.map func

    let lowerHull, upperHull = arsComputeHulls domain S fS

    arsPlot upperHull lowerHull domain S fS func

    upperHull
    |> Array.map (fun uh -> uh.Right - uh.Left, uh.Pr)
    |> Array.iter (fun (width, prob) -> printfn "%A : %A" width prob)


    let rec getSamples lowerHull upperHull S (samples: float list) =
        if samples.Length = nSamples then 
            samples
        else
            // sample x from Hull
            let x = arsSampleUpperHull upperHull rnd
            let lhVal, uhVal = arsEvalHulls x lowerHull upperHull 

            let U = rnd.NextDouble()
           
            // three cases for acception/rejection
            // computation differs from the matlab one - it was incorrect
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


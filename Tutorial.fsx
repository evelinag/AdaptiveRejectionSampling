#nowarn "211"
#r "packages/MathNet.Numerics.3.5.0/lib/net40/MathNet.Numerics.dll"
#r "packages/MathNet.Numerics.FSharp.3.5.0/lib/net40/MathNet.Numerics.FSharp.dll"
#I "packages/RProvider.1.1.6"
#load "RProvider.fsx"

open System
open RDotNet
open RProvider
open RProvider.graphics
open RProvider.grDevices

open System
open MathNet.Numerics

#load "ars.fs"
open AdaptiveRejectionSampling

//=================================================
// Sample 1-D Gaussian

// un-normalized Gaussian log likelihood
let gaussianFunc x = 
    let sigma = 1.0
    -(x*x) / sigma

let gaussianA = - 1.0
let gaussianB = 1.0
let gaussianDomain = (-infinity, infinity)

// draw samples
let samplesGaussian = 
    adaptiveRejectionSampling gaussianFunc gaussianA gaussianB gaussianDomain 10000

R.x11()
namedParams [ "x", box samplesGaussian; "breaks", box 20]
|> R.hist

// qq-plot
open RProvider.stats
R.qqnorm(y=samplesGaussian)
R.qqline(samplesGaussian)

// =============================================
// Sample Beta distribution

// un-normalized beta log likelihood
let betaFunc alpha beta x = 
    (alpha - 1.0) * log x + (beta - 1.0) * log (1.0 - x)

let alpha, beta = 2.0, 5.0
let betaA = 0.2
let betaB = 0.8
let betaDomain = (0.0, 1.0) // finite bounds
let samplesBeta = 
    adaptiveRejectionSampling (betaFunc alpha beta) betaA betaB betaDomain 10000

R.x11()
namedParams [ "x", box samplesBeta; "breaks", box 50; "probability", box true]
|> R.hist

// Compare with samples from the true distribution
open RProvider.stats
R.qqplot(samplesBeta, R.rbeta(1000, alpha, beta))

// ====================================================
// Log likelihood from a more complex model

let compositeWeights = [|0.25; 0.25; 0.25; 0.25|] 
let pis = [|0.25; 0.25; 0.25; 0.25|]
let a0 = 0.01
let b0 = 20.0

let func x = 
    (a0 - 1.0) * log x - b0 * x
    + ((Array.sum compositeWeights) * x |> SpecialFunctions.GammaLn)
    - (Array.sumBy (fun p -> SpecialFunctions.GammaLn(p * x)) compositeWeights)
    + (Array.map2 (fun prior_p current_p -> 
        (prior_p  * x - 1.0) * log current_p) compositeWeights pis
        |> Array.sum)

// Plot the function
let xs = [| 0.001 .. 0.001 .. 1.0|]
let ys = xs |> Array.map func
let ys' = ys |> Array.map exp

R.x11()
// concave log likelihood
R.plot(namedParams["x", box xs; "y", box ys; "type", box "l"])
// actual likelihood
R.plot(namedParams["x", box xs; "y", box ys'; "type", box "l"])


let domain = (0.0, infinity)
let a = 0.001
let b = 0.5
let nSamples = 10000
let samples = adaptiveRejectionSampling func a b domain nSamples

R.x11()
namedParams [ "x", box samples; "breaks", box 100; "xlim", box [0.0; 1.0]]
|> R.hist


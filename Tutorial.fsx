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

//#load "ars-rewrite.fs"
//open AdaptiveRejectionSampling.ArsRewrite

#load "ars.fs"
open AdaptiveRejectionSampling
// =======================

let compositeWeights = [|0.25; 0.25; 0.25; 0.25|] 
let pis = [|0.25; 0.25; 0.25; 0.25|]
let a0 = 2.0
let b0 = 20.0

let func x = 
    (a0 - 1.0) * log x - b0 * x
    + ((Array.sum compositeWeights) * x |> SpecialFunctions.GammaLn)
    - (Array.sumBy (fun p -> SpecialFunctions.GammaLn(p * x)) compositeWeights)
    + (Array.map2 (fun prior_p current_p -> 
        (prior_p  * x - 1.0) * log current_p) compositeWeights pis
        |> Array.sum)

let xs = [| 0.001 .. 0.001 .. 1.0|]
let ys = xs |> Array.map func
let ys' = ys |> Array.map exp

//R.x11()
R.plot(namedParams["x", box xs; "y", box ys; "type", box "l"])
R.plot(namedParams["x", box xs; "y", box ys'; "type", box "l"])

let domain = (0.0, Double.PositiveInfinity)
let a = 0.1
let b = 0.5
let nSamples = 10


let samples = adaptiveRejectionSampling func a b domain 50000

R.x11()
namedParams [ "x", box samples; "breaks", box 100; "xlim", box [0.0; 1.0]]
|> R.hist

// look at accepts-rejects
acceptsRejects |> List.sumBy (fun (x,b) -> if b = false then 1 else 0)
acceptsRejects |> List.choose (fun (x,b) -> if b = false then Some x else None) |> R.plot



let func = (fun (x:float) -> - x*x)
let a = -1.0
let b = 1.0
let domain = (Double.NegativeInfinity, Double.PositiveInfinity)
let nSamples = 1000
let quad = ars func a b domain nSamples

R.hist(quad)

let sampleAvg = 
    [| for i in 0..100 ->
        let quad = ars func a b domain nSamples
        quad |> List.average |]
    |> Array.average


//  ===============================================
// Sample 1-D Gaussian
let gaussianFunc x = 
    let sigma = 1.0
    -(x*x) / sigma

let a = - 1.0
let b = 1.0
let domain = (Double.NegativeInfinity, Double.PositiveInfinity)
let nSamples = 1000

// draw samples
let samplesGaussian = adaptiveRejectionSampling gaussianFunc a b domain nSamples

R.x11()
namedParams [ "x", box samplesGaussian; "breaks", box 20]
|> R.hist

open RProvider.StatDA
R.ppplot_das(samplesGaussian)

// =============================================
// Sample Beta distribution
let betaFunc alpha beta x = 
    (alpha - 1.0) * log x + (beta - 1.0) * log (1.0 - x)

R.x11()
[ 0.0 .. 0.1 .. 1.0] |> List.map (betaFunc 2.0 5.0) |> R.plot

let alpha, beta = 2.0, 5.0
let a = 0.2
let b = 0.8
let domain = (0.0, 1.0)
let func = betaFunc alpha beta
let samplesBeta = adaptiveRejectionSampling func a b domain 10000

R.x11()
namedParams [ "x", box samplesBeta; "breaks", box 50; "probability", box true]
|> R.hist

open RProvider.stats
R.qqplot(samplesBeta, R.rbeta(1000, alpha, beta))
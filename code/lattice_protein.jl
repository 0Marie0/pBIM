# General imports
using Revise
using StaticArrays
using DataStructures
using CSV
using DataFrames
using Random
using ProgressLogging
using ProgressMeter
using Plots
using Statistics


println("Imports done, importing hamiltonian path")
includet("file_hamiltonian_path.jl")

println("hamiltonian path imported, importing metropolis")
includet("file_metropolis.jl")

println("metropolis imported, importing benchmark")
includet("file_benchmark.jl")

println("benchmark imported, importing visualization")
includet("visualization.jl")
# Nu の分母（純拡散解の浮力分散散逸）を数値的に求めるための共通関数
#
# issue #5 の定義：Nu = ⟨χ⟩ / ⟨χ_diff⟩．χ = κ(∂x b² + ∂z b²) を公式例と同じ Oceananigans の演算子
# （@at (Center, Center, Center) と ∫dV）で求め，分母も純拡散解 b_diff に「同じ演算子」を適用した数値を使う．
# b_diff は run_horizontal_convection.jl と同じ上面の値境界条件を持つ Field に set! し，
# fill_halo_regions! で halo を埋めてから微分する（上面の勾配は halo = 2 b_s - b_N から決まる）．
#
# plot_horizontal_convection.jl と check_diffusive_limit.jl から include して使う．

using Oceananigans
using Oceananigans.BoundaryConditions: regularize_field_boundary_conditions, fill_halo_regions!

# 上面の浮力（run_horizontal_convection.jl と同じ関数・同じ引数）
@inline bˢ(x, t, p) = - p.b★ * cos(2π * x / p.Lx)

"""run スクリプトと同じ境界条件を，モデルを作らずに Field 単体で使える形にして返す．"""
function buoyancy_boundary_conditions(grid; b★, Lx)
    bcs = FieldBoundaryConditions(top = ValueBoundaryCondition(bˢ, parameters=(; b★, Lx)))
    return regularize_field_boundary_conditions(bcs, grid, :b)
end

"""純拡散解 b_diff(x,z) = b_s(x) cosh(2π(H+z)/Lx) / cosh(2πH/Lx) を halo つきで返す．"""
function diffusive_buoyancy(grid; b★, Lx, H)
    b = CenterField(grid; boundary_conditions = buoyancy_boundary_conditions(grid; b★, Lx))
    set!(b, (x, z) -> -b★ * cos(2π * x / Lx) * cosh(2π * (H + z) / Lx) / cosh(2π * H / Lx))
    fill_halo_regions!(b, Clock(time = 0.0), (;))
    return b
end

"""公式例と同じ演算子での領域平均 ⟨χ⟩ = ∫ κ(∂x b² + ∂z b²) dV / (Lx H)．"""
function mean_dissipation(b, κ; Lx, H)
    χ = @at (Center, Center, Center) κ * (∂x(b)^2 + ∂z(b)^2)
    I = Field(∫dV(χ / (Lx * H)))
    compute!(I)
    return I[1, 1, 1]
end

"""解析値 χ_diff = κ b★² π tanh(2πH/Lx) / (Lx H)（公式例の式）．"""
analytic_chi_diff(κ; b★, Lx, H) = κ * b★^2 * π * tanh(2π * H / Lx) / (Lx * H)

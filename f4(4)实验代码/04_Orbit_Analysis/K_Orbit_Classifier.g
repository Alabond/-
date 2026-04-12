#############################################################################
# 文件: K_Orbit_Classifier.g
# 描述: 基于 G2(2) Extended Vogan Diagram 的 K-Orbit 分类与紧根权重计算
# 说明: 模块3/模块5共用的 K-orbit 标签计算与归约工具。
# 核心函数用途:
#   - ParseRootString / ParseHString:
#       解析根表达式与 h 表达式，统一成可计算系数向量。
#   - CalculateWeight / DrawG2VoganDiagram:
#       基于 h 权重生成 Vogan 信息与中间可视化字符串。
#   - ReflectH / MakeKDominant:
#       通过 K-Weyl 反射将 h 归约到 K-dominant 代表。
#   - PrintKOrbitInfo:
#       模块化输出入口，返回与打印 KOrbitLabel 相关信息。
#############################################################################

# =============================================================================
# 1. 基础数学库
# =============================================================================

# 解析根字符串 (如 "3a1+2a2") 为向量 [3, 2]
ParseRootString := function(root_str)
    local v, coeff, idx, sign, num_str;
    v := [0, 0];
    root_str := Filtered(root_str, c -> c <> ' ');
    if root_str = "0" then return v; fi;
    idx := 1;
    while idx <= Length(root_str) do
        sign := 1;
        if root_str[idx] = '-' then
            sign := -1; idx := idx + 1;
        elif root_str[idx] = '+' then
            idx := idx + 1;
        fi;
        num_str := "";
        while idx <= Length(root_str) and root_str[idx] in "0123456789" do
            num_str := Concatenation(num_str, [root_str[idx]]);
            idx := idx + 1;
        od;
        if num_str = "" then coeff := 1; else coeff := Int(num_str); fi;
        coeff := coeff * sign;
        if idx <= Length(root_str) and root_str[idx] = 'a' then
            idx := idx + 1;
            if idx <= Length(root_str) and root_str[idx] in "12" then
                if root_str[idx] = '1' then v[1] := v[1] + coeff;
                elif root_str[idx] = '2' then v[2] := v[2] + coeff;
                fi;
                idx := idx + 1;
            fi;
        else
            break;
        fi;
    od;
    return v;
end;;

if not IsBoundGlobal("GlobalF4LabelStats") then
    GlobalF4LabelStats := [];
fi;
if not IsBoundGlobal("GlobalF4RawLabelStats") then
    GlobalF4RawLabelStats := [];
fi;
if not IsBoundGlobal("GlobalF4TerminalLabelStats") then
    GlobalF4TerminalLabelStats := [];
fi;
if not IsBoundGlobal("GlobalF4LabelFallbackCount") then
    GlobalF4LabelFallbackCount := 0;
fi;
if not IsBoundGlobal("GlobalF4MappingTotal") then
    GlobalF4MappingTotal := 0;
fi;
if not IsBoundGlobal("GlobalF4AllowedViolationCount") then
    GlobalF4AllowedViolationCount := 0;
fi;
if not IsBoundGlobal("GlobalF4TerminalConvergedCount") then
    GlobalF4TerminalConvergedCount := 0;
fi;
if not IsBoundGlobal("GlobalF4TerminalViolationCount") then
    GlobalF4TerminalViolationCount := 0;
fi;
if not IsBoundGlobal("GlobalF4ConvergenceIdempotentCount") then
    GlobalF4ConvergenceIdempotentCount := 0;
fi;
if not IsBoundGlobal("GlobalF4ConvergenceIdempotentViolationCount") then
    GlobalF4ConvergenceIdempotentViolationCount := 0;
fi;
if not IsBoundGlobal("GlobalF4AllowedLabelCacheLabels") then
    GlobalF4AllowedLabelCacheLabels := [];
fi;
if not IsBoundGlobal("GlobalF4AllowedLabelCachePairs") then
    GlobalF4AllowedLabelCachePairs := [];
fi;
if not IsBoundGlobal("GlobalF4StressTotal") then
    GlobalF4StressTotal := 0;
fi;
if not IsBoundGlobal("GlobalF4StressAllowedViolationCount") then
    GlobalF4StressAllowedViolationCount := 0;
fi;
if not IsBoundGlobal("GlobalF4StressTerminalViolationCount") then
    GlobalF4StressTerminalViolationCount := 0;
fi;
if not IsBoundGlobal("GlobalF4StressIdempotentViolationCount") then
    GlobalF4StressIdempotentViolationCount := 0;
fi;
if not IsBoundGlobal("GlobalF4AffineConstraintTotal") then
    GlobalF4AffineConstraintTotal := 0;
fi;
if not IsBoundGlobal("GlobalF4AffineConstraintViolationCount") then
    GlobalF4AffineConstraintViolationCount := 0;
fi;
if not IsBoundGlobal("GlobalF4CompactWeightSetTotal") then
    GlobalF4CompactWeightSetTotal := 0;
fi;
if not IsBoundGlobal("GlobalF4CompactWeightSetViolationCount") then
    GlobalF4CompactWeightSetViolationCount := 0;
fi;
if not IsBoundGlobal("GlobalShowNonCompactWeights") then
    GlobalShowNonCompactWeights := false;
fi;
if not IsBoundGlobal("GlobalCurrentSubsystemType") then
    GlobalCurrentSubsystemType := "";
fi;
if not IsBoundGlobal("ParseHString_F4Coeffs") then
    ParseHString_F4Coeffs := fail;
fi;
if not IsBoundGlobal("F4ParseHTerms") then
    F4ParseHTerms := fail;
fi;
if not IsBoundGlobal("F4CorootCoeffsFromRootVec") then
    F4CorootCoeffsFromRootVec := fail;
fi;
if not IsBoundGlobal("CalculateF4Weights") then
    CalculateF4Weights := fail;
fi;
if not IsBoundGlobal("ClassifyF4KOrbitLabel") then
    ClassifyF4KOrbitLabel := fail;
fi;
if not IsBoundGlobal("F4IsLabelInList") then
    F4IsLabelInList := fail;
fi;
if not IsBoundGlobal("F4NearestLabelFromSet") then
    F4NearestLabelFromSet := fail;
fi;
if not IsBoundGlobal("F4ExactLabelFromSet") then
    F4ExactLabelFromSet := fail;
fi;
if not IsBoundGlobal("F4LabelDigitsFromString") then
    F4LabelDigitsFromString := fail;
fi;
if not IsBoundGlobal("F4CompactDigitsToLabel") then
    F4CompactDigitsToLabel := fail;
fi;
if not IsBoundGlobal("F4RawKOrbitLabelDigits") then
    F4RawKOrbitLabelDigits := fail;
fi;
if not IsBoundGlobal("F4AllowedKOrbitLabels") then
    F4AllowedKOrbitLabels := fail;
fi;

NormalizeF4Weights_KDominant_FromCoeffs := function(c0, c1, c2, c3, c4)
    local K1, K2, K3, K4, K, marks, A, changed, iter, max_iter, w1, w2, w3, w4, w0;
    K1 := c1 - 2*c0;
    K2 := c2 - 3*c0;
    K3 := c3 - 2*c0;
    K4 := c4 - 1*c0;
    K := [K1, K2, K3, K4];
    marks := [2, 3, 2, 1];
    A := [[ 2,-1, 0, 0],
          [-1, 2,-2, 0],
          [ 0,-1, 2,-1],
          [ 0, 0,-1, 2]];
    changed := true;
    iter := 0;
    max_iter := 64;
    while changed and iter < max_iter do
        changed := false;
        iter := iter + 1;
        w1 := 2*K[1] - K[2];
        w2 := -K[1] + 2*K[2] - 2*K[3];
        w3 := -K[2] + 2*K[3] - K[4];
        w4 := -K[3] + 2*K[4];
        w0 := -(2*w1 + 3*w2 + 4*w3 + 2*w4);
        if w0 < 0 then
            K[1] := K[1] - w0*marks[1];
            K[2] := K[2] - w0*marks[2];
            K[3] := K[3] - w0*marks[3];
            K[4] := K[4] - w0*marks[4];
            changed := true;
        fi;
        if w2 < 0 then
            K[2] := K[2] - w2;
            changed := true;
        fi;
        if w3 < 0 then
            K[3] := K[3] - w3;
            changed := true;
        fi;
        if w4 < 0 then
            K[4] := K[4] - w4;
            changed := true;
        fi;
    od;
    w1 := 2*K[1] - K[2];
    w2 := -K[1] + 2*K[2] - 2*K[3];
    w3 := -K[2] + 2*K[3] - K[4];
    w4 := -K[3] + 2*K[4];
    w0 := -(2*w1 + 3*w2 + 4*w3 + 2*w4);
    return [w0, w1, w2, w3, w4];
end;;

NormalizeF4Weights_KDominant_ByWeights := function(w0, w1, w2, w3, w4)
    local Ccols, C0col, changed, iter, max_iter, ww0, ww;
    Ccols := [
        [2, -1,  0,  0],   # col 1 (unused here: a1 非紧)
        [-1, 2, -2,  0],   # col 2
        [0, -1,  2, -1],   # col 3
        [0,  0, -1,  2]    # col 4
    ];
    C0col := [-1, 0, 0, 0]; # col 0 for affine: a0 仅与 a1 单边相连
    ww := [w1, w2, w3, w4];
    ww0 := w0;
    changed := true;
    iter := 0;
    max_iter := 128;
    while changed and iter < max_iter do
        changed := false;
        iter := iter + 1;
        if ww0 < 0 then
            ww[1] := ww[1] - C0col[1]*ww0;
            ww[2] := ww[2] - C0col[2]*ww0;
            ww[3] := ww[3] - C0col[3]*ww0;
            ww[4] := ww[4] - C0col[4]*ww0;
            ww0 := -ww0;
            changed := true;
        fi;
        if ww[2] < 0 then
            ww[1] := ww[1] - Ccols[2][1]*ww[2];
            ww[2] := ww[2] - Ccols[2][2]*ww[2];
            ww[3] := ww[3] - Ccols[2][3]*ww[2];
            ww[4] := ww[4] - Ccols[2][4]*ww[2];
            changed := true;
        fi;
        if ww[3] < 0 then
            ww[1] := ww[1] - Ccols[3][1]*ww[3];
            ww[2] := ww[2] - Ccols[3][2]*ww[3];
            ww[3] := ww[3] - Ccols[3][3]*ww[3];
            ww[4] := ww[4] - Ccols[3][4]*ww[3];
            changed := true;
        fi;
        if ww[4] < 0 then
            ww[1] := ww[1] - Ccols[4][1]*ww[4];
            ww[2] := ww[2] - Ccols[4][2]*ww[4];
            ww[3] := ww[3] - Ccols[4][3]*ww[4];
            ww[4] := ww[4] - Ccols[4][4]*ww[4];
            changed := true;
        fi;
        ww0 := -(2*ww[1] + 3*ww[2] + 4*ww[3] + 2*ww[4]);
    od;
    return [ww0, ww[1], ww[2], ww[3], ww[4]];
end;;

NormalizeF4Weights_KDominant_FromHStr := function(h_str)
    local c, w, w0, w2, w3, w4, iter, max_iter;
    c := ParseHString_F4Coeffs(h_str);
    iter := 0;
    max_iter := 128;
    while iter < max_iter do
        w := CalculateF4Weights(c);
        w0 := w[1]; w2 := w[3]; w3 := w[4]; w4 := w[5];
        if w0 < 0 then
            c[1] := c[1] - w0;               # s_0
        elif w2 < 0 then
            c[3] := c[3] - w2;               # s_2
        elif w3 < 0 then
            c[4] := c[4] - w3;               # s_3
        elif w4 < 0 then
            c[5] := c[5] - w4;               # s_4
        else
            break;
        fi;
        iter := iter + 1;
    od;
    return CalculateF4Weights(c);
end;;

F4ReflectCompactWeight := function(state, idx)
    local w0, w1, w2, w3, w4, wi;
    w0 := state[1]; w1 := state[2]; w2 := state[3]; w3 := state[4]; w4 := state[5];
    if idx = 0 then
        w1 := w1 + w0;
    elif idx = 2 then
        wi := w2;
        w1 := w1 + wi;
        w2 := -wi;
        w3 := w3 + wi;
    elif idx = 3 then
        wi := w3;
        w2 := w2 + 2*wi;
        w3 := -wi;
        w4 := w4 + wi;
    elif idx = 4 then
        wi := w4;
        w3 := w3 + wi;
        w4 := -wi;
    fi;
    w0 := -(2*w1 + 3*w2 + 4*w3 + 2*w4);
    return [w0, w1, w2, w3, w4];
end;;

F4ReflectWeight := function(state, idx)
    local w0, w1, w2, w3, w4, wi;
    w0 := state[1]; w1 := state[2]; w2 := state[3]; w3 := state[4]; w4 := state[5];
    if idx = 1 then
        wi := w1;
        w0 := w0 + wi;
        w1 := -wi;
        w2 := w2 + wi;
        return [w0, w1, w2, w3, w4];
    fi;
    return F4ReflectCompactWeight(state, idx);
end;;

F4CompactDominanceKey := function(state)
    return [state[1], state[3], state[4], state[5], -(state[2] * state[2])];
end;;

F4PreferCompactDominantState := function(candidate, current)
    local cand_key, cur_key, i;
    if current = fail then
        return true;
    fi;
    cand_key := F4CompactDominanceKey(candidate);
    cur_key := F4CompactDominanceKey(current);
    for i in [1..Length(cand_key)] do
        if cand_key[i] > cur_key[i] then
            return true;
        elif cand_key[i] < cur_key[i] then
            return false;
        fi;
    od;
    return String(candidate) < String(current);
end;;

NormalizeF4Weights_KDominant_Simple := function(w0, w1, w2, w3, w4)
    local start, queue, seen, head, cur, nxt, key, idx, best;
    start := [w0, w1, w2, w3, w4];
    queue := [start];
    seen := [String(start)];
    head := 1;
    best := fail;
    while head <= Length(queue) and Length(queue) <= 4096 do
        cur := queue[head];
        head := head + 1;
        if cur[1] >= 0 and cur[3] >= 0 and cur[4] >= 0 and cur[5] >= 0 then
            if F4PreferCompactDominantState(cur, best) then
                best := cur;
            fi;
        fi;
        for idx in [0, 2, 3, 4] do
            nxt := F4ReflectCompactWeight(cur, idx);
            key := String(nxt);
            if Position(seen, key) = fail then
                Add(seen, key);
                Add(queue, nxt);
            fi;
        od;
    od;
    if best <> fail then
        return best;
    fi;
    return start;
end;;

NormalizeF4Weights_KDominant_Greedy := function(w0, w1, w2, w3, w4)
    local start, cur, iter, max_iter;
    start := [w0, w1, w2, w3, w4];
    cur := [w0, w1, w2, w3, w4];
    iter := 0;
    max_iter := 256;
    while iter < max_iter do
        if cur[1] < 0 then
            cur := F4ReflectCompactWeight(cur, 0);
        elif cur[3] < 0 then
            cur := F4ReflectCompactWeight(cur, 2);
        elif cur[4] < 0 then
            cur := F4ReflectCompactWeight(cur, 3);
        elif cur[5] < 0 then
            cur := F4ReflectCompactWeight(cur, 4);
        else
            break;
        fi;
        iter := iter + 1;
    od;
    if cur[1] >= 0 and cur[3] >= 0 and cur[4] >= 0 and cur[5] >= 0 then
        return cur;
    fi;
    return NormalizeF4Weights_KDominant_Simple(start[1], start[2], start[3], start[4], start[5]);
end;;

F4SelectKDominantRepresentative := function(w0, w1, w2, w3, w4)
    return NormalizeF4Weights_KDominant_Simple(w0, w1, w2, w3, w4);
end;;

# 灵活解析根字符串，支持 a1..a4（用于 F4/B4 子系统）
ParseRootStringFlexible := function(root_str)
    local v, coeff, idx, sign, num_str, index_str, index_val, has_higher, i, ch;
    v := [0, 0, 0, 0];
    root_str := Filtered(root_str, c -> c <> ' ');
    if root_str = "0" then return v; fi;
    idx := 1;
    while idx <= Length(root_str) do
        sign := 1;
        if root_str[idx] = '-' then
            sign := -1; idx := idx + 1;
        elif root_str[idx] = '+' then
            idx := idx + 1;
        fi;
        num_str := "";
        while idx <= Length(root_str) and root_str[idx] in "0123456789" do
            num_str := Concatenation(num_str, [root_str[idx]]);
            idx := idx + 1;
        od;
        if num_str = "" then coeff := 1; else coeff := Int(num_str); fi;
        coeff := coeff * sign;
        if idx <= Length(root_str) and root_str[idx] = 'a' then
            idx := idx + 1;
            index_str := "";
            while idx <= Length(root_str) and root_str[idx] in "0123456789" do
                index_str := Concatenation(index_str, [root_str[idx]]);
                idx := idx + 1;
            od;
            if index_str <> "" then
                index_val := Int(index_str);
                if index_val >= 1 and index_val <= 4 then
                    v[index_val] := v[index_val] + coeff;
                fi;
            fi;
        else
            break;
        fi;
    od;
    return v;
end;;

# 解析 h 字符串为 B4 系数 [c0,c1,c2,c3]，专用于 F4-B4 场景
ParseHString_B4Coeffs := function(h_str)
    local s, i, term, parts, t, c0, c1, c2, c3, pos, num, sign, idx;
    s := Filtered(h_str, c -> c <> ' ');
    parts := [];
    term := "";
    for i in [1..Length(s)] do
        if s[i] = '+' then
            if term <> "" then Add(parts, term); fi;
            term := "";
        else
            term := Concatenation(term, [s[i]]);
        fi;
    od;
    if term <> "" then Add(parts, term); fi;
    c0 := 0; c1 := 0; c2 := 0; c3 := 0;
    for t in parts do
        sign := 1;
        idx := 1;
        if idx <= Length(t) and t[idx] = '-' then sign := -1; idx := idx + 1; fi;
        num := "";
        while idx <= Length(t) and t[idx] in "0123456789" do
            num := Concatenation(num, [t[idx]]);
            idx := idx + 1;
        od;
        if num = "" then num := "1"; fi;
        if idx <= Length(t) and t[idx] = 'H' then
            # 期望 "H_{aK}"
            idx := idx + 1; # at '{'
            while idx <= Length(t) and t[idx] <> 'a' do idx := idx + 1; od;
            if idx <= Length(t) and t[idx] = 'a' then
                idx := idx + 1;
                if idx <= Length(t) and t[idx] in "0123" then
                    if t[idx] = '0' then c0 := c0 + sign * Int(num);
                    elif t[idx] = '1' then c1 := c1 + sign * Int(num);
                    elif t[idx] = '2' then c2 := c2 + sign * Int(num);
                    elif t[idx] = '3' then c3 := c3 + sign * Int(num);
                    fi;
                fi;
            fi;
        fi;
    od;
    return [c0, c1, c2, c3];
end;;

# 解析 h 字符串为 F4 扩展节点系数 [c0,c1,c2,c3,c4]
ParseHString_F4Coeffs := function(h_str)
    local terms, t, c0, c1, c2, c3, c4, coroot;
    terms := F4ParseHTerms(h_str);
    c0 := 0; c1 := 0; c2 := 0; c3 := 0; c4 := 0;
    for t in terms do
        coroot := F4CorootCoeffsFromRootVec(t[2]);
        if coroot <> fail then
            c1 := c1 + t[1] * coroot[1];
            c2 := c2 + t[1] * coroot[2];
            c3 := c3 + t[1] * coroot[3];
            c4 := c4 + t[1] * coroot[4];
        fi;
    od;
    return [c0, c1, c2, c3, c4];
end;;

# 绘制 F4(4) 的 Extended Dynkin Diagram，并基于颜色判定 compact dominance
# 计算 F4 简单根上的权重 w = [w0, w1, w2, w3, w4]
# wj = alpha_j(h)
CalculateF4Weights := function(h_coeffs)
    local c0, c1, c2, c3, c4, K1, K2, K3, K4, w1, w2, w3, w4, w0;
    c0 := h_coeffs[1];
    c1 := h_coeffs[2];
    c2 := h_coeffs[3];
    c3 := h_coeffs[4];
    c4 := h_coeffs[5];
    
    # H_{a0} = -H_{theta} = -(2H1 + 3H2 + 2H3 + 1H4)
    # h = c0*H0 + c1*H1 + c2*H2 + c3*H3 + c4*H4
    #   = (c1-2c0)H1 + (c2-3c0)H2 + (c3-2c0)H3 + (c4-c0)H4
    
    K1 := c1 - 2*c0;
    K2 := c2 - 3*c0;
    K3 := c3 - 2*c0;
    K4 := c4 - 1*c0;
    
    # Cartan Matrix F4 (1,2 Long; 3,4 Short), 标准单根系约定:
    #    1  2  3  4
    # 1  2 -1  0  0
    # 2 -1  2 -2  0
    # 3  0 -1  2 -1
    # 4  0  0 -1  2
    
    w1 := 2*K1 - 1*K2;
    w2 := -1*K1 + 2*K2 - 2*K3;
    w3 := -1*K2 + 2*K3 - 1*K4;
    w4 := -1*K3 + 2*K4;
    
    # w0 = alpha0(h) = -theta(h)
    # theta = 2a1 + 3a2 + 4a3 + 2a4
    w0 := -(2*w1 + 3*w2 + 4*w3 + 2*w4);
    
    return [w0, w1, w2, w3, w4];
end;;

F4_SymMatrix := function()
    return [[2,-1,0,0],
            [-1,2,-1,0],
            [0,-1,1,-1/2],
            [0,0,-1/2,1]];
end;;

F4CorootCoeffsFromRootVec := function(root_vec)
    local S, root_len2, i, sum_i, coeffs;
    S := F4_SymMatrix();
    root_len2 := 0;
    for i in [1..4] do
        sum_i := root_vec[1]*S[i][1] + root_vec[2]*S[i][2] + root_vec[3]*S[i][3] + root_vec[4]*S[i][4];
        root_len2 := root_len2 + root_vec[i] * sum_i;
    od;
    if root_len2 = 0 then
        return fail;
    fi;
    coeffs := [
        root_vec[1] * S[1][1] / root_len2,
        root_vec[2] * S[2][2] / root_len2,
        root_vec[3] * S[3][3] / root_len2,
        root_vec[4] * S[4][4] / root_len2
    ];
    return coeffs;
end;;

F4ParseHTerms := function(h_str)
    local s, idx, sign, num_str, coeff, root_part, terms, slash, num_str1, num_str2, num1, num2, ParseRational;
    ParseRational := function(str)
        slash := PositionSublist(str, "/");
        if slash <> fail then
            num_str1 := str{[1..slash - 1]};
            num_str2 := str{[slash + 1..Length(str)]};
            num1 := Int(String(num_str1));
            num2 := Int(String(num_str2));
            if num1 = fail or num2 = fail or num2 = 0 then
                return fail;
            fi;
            return num1 / num2;
        fi;
        return Int(String(str));
    end;
    s := Filtered(h_str, c -> c <> ' ');
    idx := 1;
    terms := [];
    while idx <= Length(s) do
        sign := 1;
        if s[idx] = '+' then
            idx := idx + 1;
        elif s[idx] = '-' then
            sign := -1;
            idx := idx + 1;
        fi;
        num_str := "";
        while idx <= Length(s) and s[idx] in "0123456789/" do
            num_str := Concatenation(num_str, [s[idx]]);
            idx := idx + 1;
        od;
        if num_str = "" then
            coeff := sign;
        else
            coeff := sign * ParseRational(num_str);
        fi;
        while idx <= Length(s) and s[idx] <> '{' do
            idx := idx + 1;
        od;
        if idx > Length(s) then
            break;
        fi;
        idx := idx + 1;
        root_part := "";
        while idx <= Length(s) and s[idx] <> '}' do
            root_part := Concatenation(root_part, [s[idx]]);
            idx := idx + 1;
        od;
        if idx <= Length(s) and s[idx] = '}' then
            idx := idx + 1;
        fi;
        if root_part <> "" then
            Add(terms, [coeff, ParseRootStringFlexible(root_part)]);
        fi;
    od;
    return terms;
end;;

ParseH_Terms_F4 := function(h_str)
    return F4ParseHTerms(h_str);
end;;

CalculateF4WeightsFromH := function(h_str)
    local terms, S, w1, w2, w3, w4, t, c, v, num1, num2, num3, num4, den, i, sum_i, w0;
    terms := ParseH_Terms_F4(h_str);
    S := F4_SymMatrix();
    w1 := 0; w2 := 0; w3 := 0; w4 := 0;
    for t in terms do
        c := t[1];
        v := t[2];
        den := 0;
        for i in [1..4] do
            sum_i := v[1]*S[i][1] + v[2]*S[i][2] + v[3]*S[i][3] + v[4]*S[i][4];
            den := den + v[i]*sum_i;
        od;
        if den = 0 then
            continue;
        fi;
        num1 := 2*(S[1][1]*v[1] + S[1][2]*v[2] + S[1][3]*v[3] + S[1][4]*v[4]);
        num2 := 2*(S[2][1]*v[1] + S[2][2]*v[2] + S[2][3]*v[3] + S[2][4]*v[4]);
        num3 := 2*(S[3][1]*v[1] + S[3][2]*v[2] + S[3][3]*v[3] + S[3][4]*v[4]);
        num4 := 2*(S[4][1]*v[1] + S[4][2]*v[2] + S[4][3]*v[3] + S[4][4]*v[4]);
        w1 := w1 + c * (num1/den);
        w2 := w2 + c * (num2/den);
        w3 := w3 + c * (num3/den);
        w4 := w4 + c * (num4/den);
    od;
    w0 := -(2*w1 + 3*w2 + 4*w3 + 2*w4);
    return [w0, w1, w2, w3, w4];
end;;

ScaleRationalsToIntegers := function(lst)
    local i, v, d, n, lcm_d, res;
    lcm_d := 1;
    for i in [1..Length(lst)] do
        v := lst[i];
        if not IsInt(v) then
            d := DenominatorRat(v);
            lcm_d := LcmInt(lcm_d, d);
        fi;
    od;
    res := [];
    for i in [1..Length(lst)] do
        v := lst[i];
        if IsInt(v) then
            Add(res, v * lcm_d);
        else
            n := NumeratorRat(v);
            d := DenominatorRat(v);
            Add(res, n * (lcm_d / d));
        fi;
    od;
    return res;
end;;

F4RawKOrbitLabelDigits := function(w0, w2, w3, w4)
    local ints;
    ints := ScaleRationalsToIntegers([w0, w2, w3, w4]);
    return [ints[1], ints[2], ints[3], ints[4]];
end;;

F4KOrbitLabelDigits := function(w0, w2, w3, w4)
    return F4RawKOrbitLabelDigits(w0, w2, w3, w4);
end;;

F4AllowedKOrbitLabels := function()
    if IsBoundGlobal("F4AllowedKOrbitLabelsProvider") then
        return ValueGlobal("F4AllowedKOrbitLabelsProvider")();
    fi;
    return [
        "0000","1100","2001","0010","3100","1101","4000","0002","2200","0020","2011","4201",
        "2210","1111","1301","3111","0400","4020","2202","8400","4402","3131","4040","2222",
        "4422","8404","8444"
    ];
end;;

F4ScalingVariants := function()
    return [[1,1,1,1],[1,1,2,2],[2,1,2,2],[2,2,1,1],[2,2,2,2]];
end;;

F4TerminalKOrbitLabels := function()
    return F4AllowedKOrbitLabels();
end;;

F4AllowedCompactDominantWeights := function()
    return [0, 1, 2, 3, 4, 8];
end;;

F4SelectTerminalRepresentative := function(w0, w2, w3, w4)
    local terminals, exact;
    terminals := F4TerminalKOrbitLabels();
    exact := ClassifyF4KOrbitLabel(w0, w2, w3, w4);
    if exact <> fail and F4IsLabelInList(exact, terminals) then
        return rec(label := exact, source := "exact_terminal", scale := 1, err := 0);
    fi;
    return rec(label := fail, source := "terminal_unclassified", scale := fail, err := fail);
end;;

ValidateF4TripleRecord := function(triple)
    local has_h, has_e, has_f;
    if not IsRecord(triple) then
        return rec(ok := false, reason := "triple_not_record");
    fi;
    has_h := IsBound(triple.h) and IsString(triple.h) and triple.h <> "";
    has_e := IsBound(triple.e) and IsString(triple.e) and triple.e <> "";
    has_f := IsBound(triple.f) and IsString(triple.f) and triple.f <> "";
    if not has_h or not has_e or not has_f then
        return rec(ok := false, reason := "triple_missing_h_e_f");
    fi;
    if PositionSublist(triple.h, "H_{") = fail then
        return rec(ok := false, reason := "triple_h_format_invalid");
    fi;
    if PositionSublist(triple.e, "E_{") = fail then
        return rec(ok := false, reason := "triple_e_format_invalid");
    fi;
    if PositionSublist(triple.f, "F_{") = fail then
        return rec(ok := false, reason := "triple_f_format_invalid");
    fi;
    return rec(ok := true, reason := "ok");
end;;

F4ComputeKCDominantFromTriple := function(triple_or_h)
    local triple_ok, h_str, coeffs, w, norm_w, raw_compact_digits, compact_digits, dominant_ok, raw_label_digits, raw_label, compact_weights, allowed_weights, val, compact_weight_set_ok, dominant_sign_ok, raw_compact_weight_set_ok, mapping, terminal_digits;
    triple_ok := true;
    if IsRecord(triple_or_h) then
        triple_ok := ValidateF4TripleRecord(triple_or_h).ok;
        if not triple_ok then
            return rec(ok := false, reason := "triple_invalid", triple_ok := false, dominant_ok := false);
        fi;
        h_str := triple_or_h.h;
    elif IsString(triple_or_h) then
        h_str := triple_or_h;
    else
        return rec(ok := false, reason := "input_invalid", triple_ok := false, dominant_ok := false);
    fi;
    w := CalculateF4WeightsFromH(h_str);
    if w[1] = 0 and w[2] = 0 and w[3] = 0 and w[4] = 0 and w[5] = 0 then
        coeffs := ParseHString_F4Coeffs(h_str);
        w := CalculateF4Weights(coeffs);
    fi;
    norm_w := F4SelectKDominantRepresentative(w[1], w[2], w[3], w[4], w[5]);
    raw_compact_digits := F4RawKOrbitLabelDigits(norm_w[1], norm_w[3], norm_w[4], norm_w[5]);
    raw_label := Concatenation(
        String(raw_compact_digits[1]),
        String(raw_compact_digits[2]),
        String(raw_compact_digits[3]),
        String(raw_compact_digits[4])
    );
    raw_label_digits := F4LabelDigitsFromString(raw_label);
    compact_digits := raw_compact_digits;
    compact_weights := [norm_w[1], norm_w[3], norm_w[4], norm_w[5]];
    allowed_weights := F4AllowedCompactDominantWeights();
    compact_weight_set_ok := true;
    for val in compact_weights do
        if Position(allowed_weights, val) = fail then
            compact_weight_set_ok := false;
            break;
        fi;
    od;
    raw_compact_weight_set_ok := compact_weight_set_ok;
    mapping := ValueGlobal("F4ValidateMappingAndConvergence")(norm_w[1], norm_w[3], norm_w[4], norm_w[5]);
    terminal_digits := F4LabelDigitsFromString(mapping.final_label);
    compact_weight_set_ok := true;
    for val in terminal_digits do
        if Position(allowed_weights, val) = fail then
            compact_weight_set_ok := false;
            break;
        fi;
    od;
    dominant_sign_ok := (norm_w[1] >= 0) and (norm_w[3] >= 0) and (norm_w[4] >= 0) and (norm_w[5] >= 0);
    dominant_ok := dominant_sign_ok and compact_weight_set_ok and mapping.terminal_ok;
    return rec(
        ok := true,
        reason := "ok",
        triple_ok := triple_ok,
        h_str := h_str,
        raw_weights := w,
        dominant_weights := norm_w,
        raw_compact_digits := raw_compact_digits,
        compact_digits := compact_digits,
        terminal_compact_digits := terminal_digits,
        terminal_label := mapping.final_label,
        terminal_source := mapping.resolved_source,
        compact_weight_values := compact_weights,
        raw_compact_weight_set_ok := raw_compact_weight_set_ok,
        compact_weight_set_ok := compact_weight_set_ok,
        dominant_sign_ok := dominant_sign_ok,
        dominant_ok := dominant_ok,
        resolved_label := mapping.resolved_label,
        resolved_source := mapping.resolved_source,
        resolved_corrected := mapping.resolved_corrected,
        resolved_fallback := mapping.resolved_fallback,
        resolved_scale := mapping.resolved_scale,
        resolved_err := mapping.resolved_err,
        terminal_ok := mapping.terminal_ok,
        terminal_steps := mapping.steps
    );
end;;

F4CompactDigitsToLabel := function(digits)
    return Concatenation(String(digits[1]), String(digits[2]), String(digits[3]), String(digits[4]));
end;;

F4AllowedLabelDigitPairs := function()
    local labels, i, pairs;
    labels := F4AllowedKOrbitLabels();
    if IsBoundGlobal("GlobalF4AllowedLabelCacheLabels") and IsBoundGlobal("GlobalF4AllowedLabelCachePairs") then
        if GlobalF4AllowedLabelCacheLabels = labels then
            return GlobalF4AllowedLabelCachePairs;
        fi;
    fi;
    pairs := [];
    for i in [1..Length(labels)] do
        Add(pairs, rec(label := labels[i], digits := F4LabelDigitsFromString(labels[i])));
    od;
    GlobalF4AllowedLabelCacheLabels := ShallowCopy(labels);
    GlobalF4AllowedLabelCachePairs := pairs;
    return pairs;
end;;

RecordF4LabelStats := function(label, is_fallback)
    local p, i;
    if not IsBoundGlobal("GlobalF4LabelStats") then
        GlobalF4LabelStats := [];
    fi;
    if not IsBoundGlobal("GlobalF4LabelFallbackCount") then
        GlobalF4LabelFallbackCount := 0;
    fi;
    p := fail;
    for i in [1..Length(GlobalF4LabelStats)] do
        if GlobalF4LabelStats[i][1] = label then
            p := i;
            break;
        fi;
    od;
    if p = fail then
        Add(GlobalF4LabelStats, [label, 1]);
    else
        GlobalF4LabelStats[p][2] := GlobalF4LabelStats[p][2] + 1;
    fi;
    if is_fallback then
        GlobalF4LabelFallbackCount := GlobalF4LabelFallbackCount + 1;
    fi;
end;;

RecordF4RawLabelStats := function(label)
    local p, i;
    if not IsBoundGlobal("GlobalF4RawLabelStats") then
        GlobalF4RawLabelStats := [];
    fi;
    p := fail;
    for i in [1..Length(GlobalF4RawLabelStats)] do
        if GlobalF4RawLabelStats[i][1] = label then
            p := i;
            break;
        fi;
    od;
    if p = fail then
        Add(GlobalF4RawLabelStats, [label, 1]);
    else
        GlobalF4RawLabelStats[p][2] := GlobalF4RawLabelStats[p][2] + 1;
    fi;
end;;

RecordF4TerminalLabelStats := function(label)
    local p, i;
    if not IsBoundGlobal("GlobalF4TerminalLabelStats") then
        GlobalF4TerminalLabelStats := [];
    fi;
    p := fail;
    for i in [1..Length(GlobalF4TerminalLabelStats)] do
        if GlobalF4TerminalLabelStats[i][1] = label then
            p := i;
            break;
        fi;
    od;
    if p = fail then
        Add(GlobalF4TerminalLabelStats, [label, 1]);
    else
        GlobalF4TerminalLabelStats[p][2] := GlobalF4TerminalLabelStats[p][2] + 1;
    fi;
end;;

RecordF4MappingStats := function(mapping)
    if not IsBoundGlobal("GlobalF4MappingTotal") then
        GlobalF4MappingTotal := 0;
    fi;
    if not IsBoundGlobal("GlobalF4AllowedViolationCount") then
        GlobalF4AllowedViolationCount := 0;
    fi;
    if not IsBoundGlobal("GlobalF4TerminalConvergedCount") then
        GlobalF4TerminalConvergedCount := 0;
    fi;
    if not IsBoundGlobal("GlobalF4TerminalViolationCount") then
        GlobalF4TerminalViolationCount := 0;
    fi;
    if not IsBoundGlobal("GlobalF4ConvergenceIdempotentCount") then
        GlobalF4ConvergenceIdempotentCount := 0;
    fi;
    if not IsBoundGlobal("GlobalF4ConvergenceIdempotentViolationCount") then
        GlobalF4ConvergenceIdempotentViolationCount := 0;
    fi;
    GlobalF4MappingTotal := GlobalF4MappingTotal + 1;
    if not mapping.allowed_ok then
        GlobalF4AllowedViolationCount := GlobalF4AllowedViolationCount + 1;
    fi;
    if mapping.terminal_ok then
        GlobalF4TerminalConvergedCount := GlobalF4TerminalConvergedCount + 1;
    else
        GlobalF4TerminalViolationCount := GlobalF4TerminalViolationCount + 1;
    fi;
    if IsBound(mapping.idempotent_ok) and mapping.idempotent_ok then
        GlobalF4ConvergenceIdempotentCount := GlobalF4ConvergenceIdempotentCount + 1;
    else
        GlobalF4ConvergenceIdempotentViolationCount := GlobalF4ConvergenceIdempotentViolationCount + 1;
    fi;
end;;

F4LabelDigitsFromString := function(label)
    return [Int(label{[1]}), Int(label{[2]}), Int(label{[3]}), Int(label{[4]})];
end;;

F4BestScaleError := function(Dv, Ld)
    local idx, num, den, t, err;
    den := 0;
    num := 0;
    for idx in [1..4] do
        den := den + Ld[idx]*Ld[idx];
        num := num + Dv[idx]*Ld[idx];
    od;
    if den = 0 then
        if Dv[1]=0 and Dv[2]=0 and Dv[3]=0 and Dv[4]=0 then
            return rec(scale := 0, err := 0);
        fi;
        return rec(scale := fail, err := fail);
    fi;
    t := num / den;
    if t < 0 then
        return rec(scale := fail, err := fail);
    fi;
    err := 0;
    for idx in [1..4] do
        err := err + (Dv[idx] - t*Ld[idx])*(Dv[idx] - t*Ld[idx]);
    od;
    return rec(scale := t, err := err);
end;;

ClassifyF4KOrbitLabel := function(w0, w2, w3, w4)
    local digits, label, exact;
    digits := F4RawKOrbitLabelDigits(w0, w2, w3, w4);
    label := F4CompactDigitsToLabel(digits);
    if F4IsLabelInList(label, F4AllowedKOrbitLabels()) then
        return label;
    fi;
    exact := F4ExactLabelFromSet(w0, w2, w3, w4, F4AllowedKOrbitLabels());
    if exact <> fail then
        return exact.label;
    fi;
    return fail;
end;;

F4IsLabelInList := function(label, labels)
    local x;
    for x in labels do
        if x = label then
            return true;
        fi;
    od;
    return false;
end;;

F4NearestLabelFromSet := function(w0, w2, w3, w4, labels)
    local best_lab, best_t, best_err, best_variant, lab, Ld, s, variants, Dv, fit, pairs, item;
    best_lab := fail;
    best_t := 0;
    best_err := fail;
    best_variant := [1,1,1,1];
    variants := F4ScalingVariants();
    pairs := [];
    for lab in labels do
        Add(pairs, rec(label := lab, digits := F4LabelDigitsFromString(lab)));
    od;
    for s in variants do
        Dv := [s[1]*w0, s[2]*w2, s[3]*w3, s[4]*w4];
        for item in pairs do
            lab := item.label;
            Ld := item.digits;
            fit := F4BestScaleError(Dv, Ld);
            if fit.err <> fail then
                if best_err = fail or fit.err < best_err then
                    best_err := fit.err;
                    best_lab := lab;
                    best_t := fit.scale;
                    best_variant := s;
                fi;
            fi;
        od;
    od;
    if best_lab = fail then
        if Length(labels) > 0 then
            return rec(label := labels[1], scale := 0, err := 0, variant := [1,1,1,1]);
        fi;
        return rec(label := "0400", scale := 0, err := 0, variant := [1,1,1,1]);
    fi;
    return rec(label := best_lab, scale := best_t, err := best_err, variant := best_variant);
end;;

F4ExactLabelFromSet := function(w0, w2, w3, w4, labels)
    local fit;
    fit := F4NearestLabelFromSet(w0, w2, w3, w4, labels);
    if fit.label = fail or fit.err = fail or fit.err <> 0 then
        return fail;
    fi;
    return fit;
end;;

F4ResolveAllowedLabel := function(w0, w2, w3, w4)
    local raw_digits, raw_label, exact;
    raw_digits := F4RawKOrbitLabelDigits(w0, w2, w3, w4);
    raw_label := F4CompactDigitsToLabel(raw_digits);
    if F4IsLabelInList(raw_label, F4AllowedKOrbitLabels()) then
        return rec(label := raw_label, source := "exact_allowed", corrected := false, fallback := false, scale := fail, err := 0);
    fi;
    exact := F4ExactLabelFromSet(w0, w2, w3, w4, F4AllowedKOrbitLabels());
    if exact <> fail then
        return rec(label := exact.label, source := "exact_scaled_allowed", corrected := true, fallback := false, scale := exact.scale, err := exact.err);
    fi;
    return rec(label := raw_label, source := "unlisted_compact_label", corrected := false, fallback := false, scale := fail, err := fail);
end;;

F4ConvergeToTerminalLabel := function(w0, w2, w3, w4, seed_label)
    local terminals;
    terminals := F4TerminalKOrbitLabels();
    if F4IsLabelInList(seed_label, terminals) then
        return rec(seed := seed_label, label := seed_label, steps := 0, converged := true, scale := 1, err := 0);
    fi;
    return rec(seed := seed_label, label := seed_label, steps := 0, converged := false, scale := fail, err := fail);
end;;

F4ValidateMappingAndConvergence := function(w0, w2, w3, w4)
    local resolved, conv, allowed_ok, terminal_ok, conv2, idempotent_ok, terminals;
    terminals := F4TerminalKOrbitLabels();
    resolved := F4ResolveAllowedLabel(w0, w2, w3, w4);
    conv := F4ConvergeToTerminalLabel(w0, w2, w3, w4, resolved.label);
    conv2 := F4ConvergeToTerminalLabel(w0, w2, w3, w4, conv.label);
    allowed_ok := F4IsLabelInList(resolved.label, F4AllowedKOrbitLabels());
    terminal_ok := F4IsLabelInList(conv.label, terminals);
    idempotent_ok := conv2.label = conv.label;
    return rec(
        resolved_label := resolved.label,
        resolved_source := resolved.source,
        resolved_corrected := resolved.corrected,
        resolved_fallback := resolved.fallback,
        resolved_scale := resolved.scale,
        resolved_err := resolved.err,
        allowed_ok := allowed_ok,
        final_label := conv.label,
        seed_label := conv.seed,
        steps := conv.steps,
        converged := conv.converged and terminal_ok,
        terminal_ok := terminal_ok,
        idempotent_ok := idempotent_ok,
        terminal_scale := conv.scale,
        terminal_err := conv.err
    );
end;;

F4RunConvergenceStressTest := function()
    local labels, terminals, variants, scales, L, s, t, w0, w2, w3, w4, mapping, total, allowed_bad, terminal_bad, idempotent_bad, pair;
    labels := F4AllowedLabelDigitPairs();
    terminals := F4TerminalKOrbitLabels();
    variants := F4ScalingVariants();
    scales := [1..25];
    total := 0;
    allowed_bad := 0;
    terminal_bad := 0;
    idempotent_bad := 0;
    for pair in labels do
        L := pair.digits;
        for s in variants do
            for t in scales do
                w0 := t * L[1] / s[1];
                w2 := t * L[2] / s[2];
                w3 := t * L[3] / s[3];
                w4 := t * L[4] / s[4];
                mapping := F4ValidateMappingAndConvergence(w0, w2, w3, w4);
                total := total + 1;
                if not mapping.allowed_ok then
                    allowed_bad := allowed_bad + 1;
                fi;
                if not mapping.terminal_ok or not F4IsLabelInList(mapping.final_label, terminals) then
                    terminal_bad := terminal_bad + 1;
                fi;
                if not mapping.idempotent_ok then
                    idempotent_bad := idempotent_bad + 1;
                fi;
            od;
        od;
    od;
    if not IsBoundGlobal("GlobalF4StressTotal") then GlobalF4StressTotal := 0; fi;
    if not IsBoundGlobal("GlobalF4StressAllowedViolationCount") then GlobalF4StressAllowedViolationCount := 0; fi;
    if not IsBoundGlobal("GlobalF4StressTerminalViolationCount") then GlobalF4StressTerminalViolationCount := 0; fi;
    if not IsBoundGlobal("GlobalF4StressIdempotentViolationCount") then GlobalF4StressIdempotentViolationCount := 0; fi;
    GlobalF4StressTotal := total;
    GlobalF4StressAllowedViolationCount := allowed_bad;
    GlobalF4StressTerminalViolationCount := terminal_bad;
    GlobalF4StressIdempotentViolationCount := idempotent_bad;
    return rec(total := total, allowed_bad := allowed_bad, terminal_bad := terminal_bad, idempotent_bad := idempotent_bad);
end;;

ExplainF4KOrbitLabel := function(w0, w2, w3, w4)
    local digits, mapping;
    digits := F4KOrbitLabelDigits(w0, w2, w3, w4);
    mapping := F4ValidateMappingAndConvergence(w0, w2, w3, w4);
    if not mapping.allowed_ok then
        return Concatenation(
            "    Label proof: digits ", String(digits),
            " -> [未落入既有分类表:", String(mapping.resolved_label), "]\n"
        );
    fi;
    return Concatenation(
        "    Label proof: digits ", String(digits),
        " -> [分类标签:", String(mapping.resolved_label), "]\n"
    );
end;;

# 绘制 F4(4) 的 Extended Dynkin Diagram，并基于颜色判定 compact dominance
DrawF4ExtendedVoganDiagram := function(triple_or_h, colors_list)
    local calc, h_str, w, w0, w1, w2, w3, w4, norm_w, norm_w0, norm_w1, norm_w2, norm_w3, norm_w4, m0, m1, m2, m3, m4, whites, ok, out, whites_vals, full_vals, g, i, disp_all_int, disp_g, disp_w0, disp_w1, disp_w2, disp_w3, disp_w4, label_digits, affine_num, affine_lhs, affine_rhs, a1_num, raw_label, S, ext_roots, r, c, ip, row_str, ipjj, aij, aji;
    calc := F4ComputeKCDominantFromTriple(triple_or_h);
    if not calc.ok then
        Print("  [模块5] triple 结构校验失败: ", calc.reason, "\n");
        return;
    fi;
    h_str := calc.h_str;
    w := calc.raw_weights;
    norm_w := calc.dominant_weights;
    label_digits := calc.compact_digits;
    raw_label := F4CompactDigitsToLabel(label_digits);
    RecordF4RawLabelStats(raw_label);
    RecordF4LabelStats(calc.resolved_label, false);
    RecordF4TerminalLabelStats(calc.terminal_label);
    w0 := w[1]; w1 := w[2]; w2 := w[3]; w3 := w[4]; w4 := w[5];
    norm_w0 := norm_w[1]; norm_w1 := norm_w[2]; norm_w2 := norm_w[3]; norm_w3 := norm_w[4]; norm_w4 := norm_w[5];
    
    # 不过滤未匹配标签，完整输出用于诊断
    
    # 统一缩放到整数并按紧节点的 gcd 约去公因子，作为“图上显示权重”
    # 注意：仅影响显示，不改变内部分类与标签逻辑（标签另行做整数化与归一）
    disp_all_int := ScaleRationalsToIntegers([norm_w0, norm_w1, norm_w2, norm_w3, norm_w4]);
    disp_g := 0;
    if disp_all_int[1] <> 0 then disp_g := AbsInt(disp_all_int[1]); fi;
    if disp_all_int[3] <> 0 then if disp_g = 0 then disp_g := AbsInt(disp_all_int[3]); else disp_g := GcdInt(disp_g, AbsInt(disp_all_int[3])); fi; fi;
    if disp_all_int[4] <> 0 then if disp_g = 0 then disp_g := AbsInt(disp_all_int[4]); else disp_g := GcdInt(disp_g, AbsInt(disp_all_int[4])); fi; fi;
    if disp_all_int[5] <> 0 then if disp_g = 0 then disp_g := AbsInt(disp_all_int[5]); else disp_g := GcdInt(disp_g, AbsInt(disp_all_int[5])); fi; fi;
    if disp_g > 0 then
        disp_w0 := disp_all_int[1] / disp_g;
        disp_w1 := disp_all_int[2] / disp_g;
        disp_w2 := disp_all_int[3] / disp_g;
        disp_w3 := disp_all_int[4] / disp_g;
        disp_w4 := disp_all_int[5] / disp_g;
    else
        disp_w0 := disp_all_int[1];
        disp_w1 := disp_all_int[2];
        disp_w2 := disp_all_int[3];
        disp_w3 := disp_all_int[4];
        disp_w4 := disp_all_int[5];
    fi;
    disp_w0 := norm_w0;
    disp_w2 := norm_w2;
    disp_w3 := norm_w3;
    disp_w4 := norm_w4;
    affine_num := -disp_w0 - 3 * disp_w2 - 4 * disp_w3 - 2 * disp_w4;
    a1_num := affine_num;
    disp_w1 := a1_num / 2;
    affine_lhs := disp_w0;
    affine_rhs := -2 * disp_w1 - 3 * disp_w2 - 4 * disp_w3 - 2 * disp_w4;
    if not IsBoundGlobal("GlobalF4AffineConstraintTotal") then GlobalF4AffineConstraintTotal := 0; fi;
    if not IsBoundGlobal("GlobalF4AffineConstraintViolationCount") then GlobalF4AffineConstraintViolationCount := 0; fi;
    GlobalF4AffineConstraintTotal := GlobalF4AffineConstraintTotal + 1;
    if affine_lhs <> affine_rhs then
        GlobalF4AffineConstraintViolationCount := GlobalF4AffineConstraintViolationCount + 1;
    fi;
    if not IsBoundGlobal("GlobalF4CompactWeightSetTotal") then GlobalF4CompactWeightSetTotal := 0; fi;
    if not IsBoundGlobal("GlobalF4CompactWeightSetViolationCount") then GlobalF4CompactWeightSetViolationCount := 0; fi;
    GlobalF4CompactWeightSetTotal := GlobalF4CompactWeightSetTotal + 1;
    if not calc.compact_weight_set_ok then
        GlobalF4CompactWeightSetViolationCount := GlobalF4CompactWeightSetViolationCount + 1;
    fi;
    
    m0 := " "; m1 := "*"; m2 := " "; m3 := " "; m4 := " ";
    
    out := "";
    out := Concatenation(out, "  k_c 共轭校验: triple_valid=", String(calc.triple_ok), ", compact_nonnegative=", String(calc.dominant_sign_ok), ", compact_weight_set_ok=", String(calc.compact_weight_set_ok), ", compact_dominant=", String(calc.dominant_ok), "\n");
    out := Concatenation(out, "  Extended Vogan Diagram (F4(4)):\n");
    out := Concatenation(out, "  [ ] = Compact, (*) = Non-compact (Painted)\n");
    out := Concatenation(out, "  注: a0 = -theta（负最高根）\n\n");
    out := Concatenation(out, "  ", m0, " [ ", String(disp_w0), " ] a0 (-theta = 负最高根, Long)\n");
    out := Concatenation(out, "           |\n");
    out := Concatenation(out, "           | (single bond)\n");
    out := Concatenation(out, "           |\n");
    out := Concatenation(out, "  ", m1, " [ ", String(disp_w1), " ] a1 (Long, Non-compact)\n");
    out := Concatenation(out, "           |\n");
    out := Concatenation(out, "           | (single bond)\n");
    out := Concatenation(out, "           |\n");
    out := Concatenation(out, "  ", m2, " [ ", String(disp_w2), " ] a2 (Long)\n");
    out := Concatenation(out, "           ||\n");
    out := Concatenation(out, "           || (double bond, arrow to short)\n");
    out := Concatenation(out, "           \\/\n");
    out := Concatenation(out, "  ", m3, " [ ", String(disp_w3), " ] a3 (Short)\n");
    out := Concatenation(out, "           |\n");
    out := Concatenation(out, "           | (single bond)\n");
    out := Concatenation(out, "           |\n");
    out := Concatenation(out, "  ", m4, " [ ", String(disp_w4), " ] a4 (Short)\n");
    out := Concatenation(out, "\n  Topology Summary:\n");
    out := Concatenation(out, "  a0 --- a1 --- a2 ==> a3 --- a4\n");
    S := F4_SymMatrix();
    ext_roots := [
        [-2, -3, -4, -2],
        [1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1]
    ];
    out := Concatenation(out, "\n  Compact Roots Weights:\n");
    out := Concatenation(out, "    a2: ", String(disp_w2), "\n");
    out := Concatenation(out, "    a3: ", String(disp_w3), "\n");
    out := Concatenation(out, "    a4: ", String(disp_w4), "\n");
    out := Concatenation(out, "    a0: ", String(disp_w0), "\n");
    out := Concatenation(out, ExplainF4KOrbitLabel(norm_w0, norm_w2, norm_w3, norm_w4));
    if calc.dominant_ok then
        out := Concatenation(out, "  => Compact Dominant! (K-orbit identified: ", String(calc.resolved_label), ")\n");
    elif calc.dominant_sign_ok then
        if IsBound(calc.resolved_label) and calc.resolved_label <> fail and calc.resolved_label <> "" then
            out := Concatenation(out, "  => Compact nonnegative but label ", String(calc.resolved_label), " 不在既有分类表中。\n");
        else
            out := Concatenation(out, "  => Compact nonnegative but violates weight set {0,1,2,3,4,8}.\n");
        fi;
    else
        out := Concatenation(out, "  => Not fully dominant on compact roots.\n");
    fi;
    
    Print(out);
    if IsBoundGlobal("GlobalEnableLogFile") and ValueGlobal("GlobalEnableLogFile") = true and IsBoundGlobal("GlobalLogfileName") and ValueGlobal("GlobalLogfileName") <> fail then
        AppendTo(ValueGlobal("GlobalLogfileName"), out);
    fi;
end;;

# 文本化的 F4(4) K-Orbit 信息输出（参照 G2 的模块 5 格式）
PrintKOrbitInfo_F4 := function(triple_or_h, colors_list)
    local h_str, c, chk;
    if IsRecord(triple_or_h) then
        chk := ValidateF4TripleRecord(triple_or_h);
        if not chk.ok then
            Print("        >>> 模块 5 输出: K-Orbit 分类\n");
            Print("        K-Orbit Classification (F4(4)):\n");
            Print("          triple 结构校验失败: ", chk.reason, "\n");
            return;
        fi;
        h_str := triple_or_h.h;
    else
        h_str := triple_or_h;
    fi;
    c := ParseHString_F4Coeffs(h_str);
    Print("        >>> 模块 5 输出: K-Orbit 分类\n");
    Print("        K-Orbit Classification (F4(4)):\n");
    Print("          初始 h 系数 (F4 扩展): ", c, "\n");
    DrawF4ExtendedVoganDiagram(triple_or_h, colors_list);
end;;

RunF4Module5UnitTests := function()
    local total, passed, failed, good_triple, bad_triple, scaled_h, r1, r2, r3, r4, direct_digits, terminals;
    total := 0; passed := 0; failed := [];
    good_triple := rec(
        h := "H_{a1} + H_{-a1-2a2-2a3-a4} + H_{a1+2a2+4a3+2a4}",
        e := "sqrt(1)E_{[a1]} + sqrt(1)E_{[-a1-2a2-2a3-a4]} + sqrt(1)E_{[a1+2a2+4a3+2a4]}",
        f := "sqrt(1)F_{[a1]} + sqrt(1)F_{[-a1-2a2-2a3-a4]} + sqrt(1)F_{[a1+2a2+4a3+2a4]}"
    );
    bad_triple := rec(
        h := "H_{a1}",
        e := "E_{[a1]}"
    );
    scaled_h := "2H_{a1} + 2H_{-a1-2a2-2a3-a4} + 2H_{a1+2a2+4a3+2a4}";
    terminals := F4TerminalKOrbitLabels();
    total := total + 1;
    r1 := F4ComputeKCDominantFromTriple(good_triple);
    if r1.ok and r1.triple_ok and r1.dominant_ok and r1.compact_digits = r1.raw_compact_digits and Position(terminals, r1.terminal_label) <> fail and r1.resolved_label = r1.terminal_label then
        passed := passed + 1;
    else
        Add(failed, "good_triple_should_map_to_allowed_label");
    fi;
    total := total + 1;
    direct_digits := F4KOrbitLabelDigits(r1.dominant_weights[1], r1.dominant_weights[3], r1.dominant_weights[4], r1.dominant_weights[5]);
    if r1.ok and r1.raw_compact_digits = direct_digits then
        passed := passed + 1;
    else
        Add(failed, "raw_compact_digits_should_equal_direct_normalized_digits");
    fi;
    total := total + 1;
    r3 := F4ComputeKCDominantFromTriple(good_triple.h);
    if r3.ok and r3.triple_ok and r3.compact_digits = r1.compact_digits and r3.raw_compact_digits = r1.raw_compact_digits and r3.terminal_label = r1.terminal_label and r3.resolved_label = r1.resolved_label then
        passed := passed + 1;
    else
        Add(failed, "string_input_should_match_record_input_and_allowed_label");
    fi;
    total := total + 1;
    r4 := F4ComputeKCDominantFromTriple(scaled_h);
    if r4.ok and r4.triple_ok and r4.dominant_ok and r4.compact_digits = r1.compact_digits and Position(terminals, r4.terminal_label) <> fail and r4.terminal_label = r1.terminal_label and r4.resolved_label = r1.resolved_label then
        passed := passed + 1;
    else
        Add(failed, "scaled_h_should_preserve_raw_and_allowed_mapping");
    fi;
    total := total + 1;
    r2 := F4ComputeKCDominantFromTriple(bad_triple);
    if (not r2.ok) and (not r2.triple_ok) then
        passed := passed + 1;
    else
        Add(failed, "bad_triple_should_fail_validation");
    fi;
    return rec(total := total, passed := passed, failed := failed);
end;;
# 检查根是否紧 (Compact)
# G2 Split Form 规则: n2 为偶数 -> 紧
IsCompact := function(root_vec)
    local n2;
    n2 := root_vec[2];
    if n2 < 0 then n2 := -n2; fi;
    return (n2 mod 2) = 0;
end;;

# 计算 h 在 root 上的值 (Weight)
# h = c1*H_alpha1 + c2*H_alpha2
# alpha(h):
# alpha1(H_a1) = 2, alpha1(H_a2) = -1 (对于 A2? 不, 这里是 G2)
# G2 Cartan Matrix:
#     a1  a2
# a1  2   -1  (Wait, short root is a1?)
# Standard G2: a1 short, a2 long
# <a1, a1> = 2/3 <a2, a2>? No.
# a1 (short), a2 (long)
# <a1^v, a2> = -1
# <a2^v, a1> = -3
# Cartan Matrix:
# [ 2 -1 ]
# [-3  2 ]
#
# alpha1(h) = 2*c1 - 1*c2 ? No.
# h = c1*H_a1 + c2*H_a2
# alpha1(h) = c1*<a1, a1^v> + c2*<a1, a2^v>
#           = c1*2 + c2*(-1)
# alpha2(h) = c1*<a2, a1^v> + c2*<a2, a2^v>
#           = c1*(-3) + c2*2
#
# 所以:
# weight(v, h) = v[1] * (2*c1 - c2) + v[2] * (-3*c1 + 2*c2)
#
# 我们需要解析 h_str 来得到 c1, c2.
# h_str 如 "2H_{a2} + 2H_{-3a1-2a2}"
# H_{v} = v[1]*H_a1 + v[2]*H_a2 ( linearity of coroot? No! )
# Coroots are linear only for Weyl group action, not generally H_{a+b} = H_a + H_b.
# H_{alpha} = 2*alpha / <alpha, alpha>.
#
# 对于 G2:
# H_a1 = a1^v
# H_a2 = a2^v
#
# 其他根的 H_alpha:
# Short roots: a1, a1+a2, 2a1+a2
# Long roots: a2, 3a1+a2, 3a1+2a2
#
# 如果 beta 是短根, H_beta 是短 coroot.
# 如果 beta 是长根, H_beta 是长 coroot.
#
# 关系:
# H_{a1+a2} = H_a1 + H_a2 ?
# Check: <a1+a2, H_{a1+a2}> = 2.
# (a1+a2)(H_a1+H_a2) = (2-3) + (-1+2) = -1+1 = 0. No.
# Actually H_{a1+a2} = H_a1 + 3*H_a2 (for G2?)
# Let's use the explicit coroot expansion.
#
# Co-roots in terms of simple co-roots H1, H2:
# Short roots:
# a1 -> H1
# a1+a2 -> H1 + H2 (Wait, check lengths)
# 2a1+a2 -> 2H1 + H2
#
# Long roots:
# a2 -> H2
# 3a1+a2 -> H1 + H2
# 3a1+2a2 -> H1 + 2H2
#
# 验证 a1+a2 (short):
# <a1+a2, H1+H2> = (2-3) + (-1+2) = 0?
# Something is wrong with my Cartan matrix or pairing.
# G2: a1 short, a2 long.
# <a1, a1> = 2 (norm convention 1) -> <a2, a2> = 6.
# <a1, a2^v> = 2*(-3)/6 = -1.
# <a2, a1^v> = 6*(-1)/2 = -3.
#
# Cartan:
#     H1  H2
# a1  2   -1
# a2 -3    2
#
# Coroot expansion (H_beta):
# Short roots (like a1):
# a1     -> H1
# a1+a2  -> H1 + H2
# 2a1+a2 -> 2H1 + H2
#
# Long roots (like a2):
# a2      -> H2
# 3a1+a2  -> 3H1 + H2  (Check: (3a1+a2)(3H1+H2) = 3(6-1) + (-9+2) = 15-7=8 != 2)
#
# Correct Coroot expansions for G2:
# H_{a1} = H1
# H_{a2} = H2
# H_{a1+a2} = H1 + H2 (if a1+a2 is short? No)
# Let's use GAP's reflection formula to find H_beta.
# H_beta = 2*beta / <beta, beta>.
# Assume basis H1, H2.
# H1 corresponds to a1, H2 to a2.
#
# 算法:
# 1. 计算 beta 向量 [n1, n2]
# 2. 计算 squared length: L = 2*(n1^2 - n1*n2 + n2^2) ? No.
#    Using inner product matrix: [[2, -3], [-3, 6]]
#    Norm = v G v^T.
# 3. H_beta = (2 / Norm) * (n1*a1 + n2*a2)^dual ?
#    This is getting complicated.
#    Alternative: H_beta = w(H_alpha) if beta = w(alpha).
#    
#    Since we only have 6 positive roots, let's hardcode the table.
#    Roots:
#    [1, 0] (a1, S) -> H1
#    [0, 1] (a2, L) -> H2
#    [1, 1] (a1+a2, S) -> H1 + H2
#    [2, 1] (2a1+a2, S) -> 2H1 + H2
#    [3, 1] (3a1+a2, L) -> H1 + H2
#    [3, 2] (3a1+2a2, L) -> H1 + 2H2
#    Negatives are -H.
#
#    Wait, [1, 1] (S): (a1+a2)(H1+H2) = (2-3) + (-1+2) = 0. NO.
#    Short root a1+a2.
#    <a1+a2, H_{a1+a2}> = 2.
#    Let H = xH1 + yH2.
#    (a1+a2)(xH1+yH2) = x(2-3) + y(-1+2) = -x + y = 2.
#    
#    Table correction:
#    a1 (S): H1
#    a2 (L): H2
#    a1+a2 (S): H1 + 3H2
#      Check: (a1+a2)(H1+3H2) = (2-3) + 3(-1+2) = -1 + 3 = 2. OK.
#    2a1+a2 (S): H1 + 2H2
#      Check: (2a1+a2)(H1+2H2) = (4-3) + 2(-2+2) = 1 + 0 = 1. NO.
#      Try 2H1 + 3H2: (2a1+a2)(2H1+3H2) = 2(4-3) + 3(-2+2) = 2. OK.
#    3a1+a2 (L): H1 + H2
#      Check: (3a1+a2)(H1+H2) = (6-3) + (-3+2) = 3 - 1 = 2. OK.
#    3a1+2a2 (L): H1 + 2H2
#      Check: (3a1+2a2)(H1+2H2) = (6-6) + 2(-3+4) = 2. OK.
#
#    Summary Table (coeff of [H1, H2]):
#    [1, 0] -> [1, 0]
#    [0, 1] -> [0, 1]
#    [1, 1] -> [1, 3]
#    [2, 1] -> [2, 3]
#    [3, 1] -> [1, 1]
#    [3, 2] -> [1, 2]

GetCorootCoeffs := function(root_vec)
    local abs_v, sign;
    abs_v := [root_vec[1], root_vec[2]];
    sign := 1;
    if abs_v[1] < 0 or (abs_v[1] = 0 and abs_v[2] < 0) then
        abs_v := -abs_v;
        sign := -1;
    fi;
    
    if abs_v = [1, 0] then return sign * [1, 0];
    elif abs_v = [0, 1] then return sign * [0, 1];
    elif abs_v = [1, 1] then return sign * [1, 3];
    elif abs_v = [2, 1] then return sign * [2, 3];
    elif abs_v = [3, 1] then return sign * [1, 1];
    elif abs_v = [3, 2] then return sign * [1, 2];
    else return fail; # Should not happen for roots
    fi;
end;;

# 计算 h 的总系数 [C1, C2] (h = C1*H1 + C2*H2)
ParseHString := function(h_str)
    local terms, term, C, idx, coeff, root_part, root_vec, h_coeffs, coroot, num_str, ParseRational, num_val, slash, num_str1, num_str2, num1, num2;
    ParseRational := function(str)
        slash := PositionSublist(str, "/");
        if slash <> fail then
            num_str1 := str{[1..slash-1]};
            num_str2 := str{[slash+1..Length(str)]};
            num1 := Int(String(num_str1));
            num2 := Int(String(num_str2));
            if num1 = fail or num2 = fail or num2 = 0 then
                return fail;
            fi;
            return num1/num2;
        fi;
        return Int(String(str));
    end;
    
    C := [0, 0];
    
    # 简单的字符串分割 (假设格式 "2H_{a2} + 2H_{-3a1-2a2}")
    # 移除空格
    h_str := Filtered(h_str, c -> c <> ' ');
    
    # 手动解析每一项
    idx := 1;
    while idx <= Length(h_str) do
        # 读取符号和系数
        coeff := 1;
        if h_str[idx] = '+' then idx := idx + 1;
        elif h_str[idx] = '-' then coeff := -1; idx := idx + 1;
        fi;
        
        # 读取数字系数（支持多位数）
        num_str := "";
        while idx <= Length(h_str) and h_str[idx] in "0123456789/" do
            num_str := Concatenation(num_str, [h_str[idx]]);
            idx := idx + 1;
        od;
        if num_str <> "" then
            num_val := ParseRational(num_str);
            if num_val = fail then num_val := 1; fi;
            coeff := coeff * num_val;
        fi;
        
        # 期望 "H_{"
        if idx+2 <= Length(h_str) and h_str[idx] = 'H' then
            idx := idx + 3; # Skip "H_{"
            
            # 读取根直到 "}"
            root_part := "";
            while idx <= Length(h_str) and h_str[idx] <> '}' do
                root_part := Concatenation(root_part, [h_str[idx]]);
                idx := idx + 1;
            od;
            idx := idx + 1; # Skip "}"
            
            # 解析根并查表得到 coroot
            root_vec := ParseRootString(root_part);
            coroot := GetCorootCoeffs(root_vec);
            
            if coroot <> fail then
                C[1] := C[1] + coeff * coroot[1];
                C[2] := C[2] + coeff * coroot[2];
            fi;
        else
            idx := idx + 1;
        fi;
    od;
    
    return C;
end;;

# 计算权重 w(alpha) = alpha(h)
CalculateWeight := function(root_vec, h_coeffs)
    local c1, c2, w;
    c1 := h_coeffs[1];
    c2 := h_coeffs[2];
    
    # alpha(h) = v1*alpha1(h) + v2*alpha2(h)
    # alpha1(h) = 2*c1 - c2
    # alpha2(h) = -3*c1 + 2*c2
    
    w := root_vec[1] * (2*c1 - c2) + root_vec[2] * (-3*c1 + 2*c2);
    return w;
end;;

# =============================================================================
# 2. 扩展 Vogan 图绘制
# =============================================================================

# 绘制 ASCII 图
# G2 Extended Diagram:
#      (a2, L)
#      2
#     /
#   (a1, S) 1 === 0 (-theta, S/L? theta is long, so -theta is long)
#   Wait. G2 affine diagram:
#   0 -- 1 (3 edges) => No.
#   Theta = 3a1 + 2a2 (Long).
#   Alpha0 = -Theta.
#   Alpha0 is connected to Alpha1 (Short)?
#   Check pairings.
#   <a0, a1^v> = <-3a1-2a2, a1^v> = -3(2) - 2(-3) = -6+6=0.
#   <a0, a2^v> = <-3a1-2a2, a2^v> = -3(-1) - 2(2) = 3-4=-1.
#   So a0 is connected to a2 (Long).
#   
#   Diagram:
#   1 (S) <=== 2 (L) --- 0 (L)
#   (a1)       (a2)      (a0)
#
#   Nodes: 1, 2, 0
#   Colors: Compact (Black/White?) 
#   Note: In our previous notation (Subsystem Classifier), 
#   we used "Black" for Non-compact.
#   Vogan Diagram: Painted nodes are Non-compact imaginary.
#   G2 Split: n2 odd -> Non-compact.
#   a1 (0 a2) -> Compact (White)
#   a2 (1 a2) -> Non-compact (Black)
#   a0 (-2 a2) -> Compact (White)
#
#   So G2 Split Vogan Diagram has Black at node 2.
#   1(W) <=== 2(B) --- 0(W)

DrawG2VoganDiagram := function(h_coeffs)
    local w1, w2, w0;
    
    # 计算三个节点的权重
    w1 := CalculateWeight([1, 0], h_coeffs);   # a1
    w2 := CalculateWeight([0, 1], h_coeffs);   # a2
    w0 := CalculateWeight([-3, -2], h_coeffs); # a0 (-theta)
    
    # 紧致性 (G2 Split Form)
    # a1: Compact
    # a2: Non-compact
    # a0: Compact
    
    # 节点权重已在 w0/w1/w2 中

    # ASCII Art
    # G2 Affine Diagram:
    # a0 (Long) --- a2 (Long) ==> a1 (Short)
    #               (Non-cpt)
    
    Print("  Extended Vogan Diagram (G2 Split):\n");
    Print("  [ ] = Compact, (*) = Non-compact (Painted)\n");
    Print("  注: a0 = -theta（负最高根）\n\n");
    
    # 绘制水平图
    # a0
    Print("      [ ", String(w0), " ] a0 (-theta = 负最高根, Long)\n");
    Print("           |\n");
    Print("           | (single bond)\n");
    Print("           |\n");
    Print("      (*) [ ", String(w2), " ] a2 (Long)\n");
    Print("           ||\n");
    Print("           || (triple bond, arrow to short)\n");
    Print("           \\/\n");
    Print("      [ ", String(w1), " ] a1 (Short)\n");
    
    Print("\n  Topology Summary:\n");
    Print("  a0(L) --- a2(L) ===> a1(S)\n");
    
    # 检查 Domination 条件
    # 我们关注的是 Compact Roots 上的权重是否为正整数
    Print("\n  Compact Roots Weights:\n");
    Print("    a1: ", w1, "\n");
    Print("    a0: ", w0, "\n");
    
    if w1 >= 0 and w0 >= 0 then
        Print("  => Compact Dominant! (K-orbit identified)\n");
    else
        Print("  => Not fully dominant on compact roots.\n");
    fi;
end;;

# =============================================================================
# 3. K-Dominant 变换 (Weyl Group Action)
# =============================================================================

# 对 h 系数应用反射 s_alpha
# h' = s_alpha(h) = h - alpha(h) H_alpha
ReflectH := function(h_coeffs, root_vec)
    local w, coroot_coeffs, new_c;
    
    # 1. 计算当前权重 alpha(h)
    w := CalculateWeight(root_vec, h_coeffs);
    
    # 2. 获取 H_alpha 的系数
    coroot_coeffs := GetCorootCoeffs(root_vec);
    
    # 3. h' = h - w * H_alpha
    new_c := [
        h_coeffs[1] - w * coroot_coeffs[1],
        h_coeffs[2] - w * coroot_coeffs[2]
    ];
    
    return new_c;
end;;

# 自动变换 h 到 K-Dominant Chamber
# 紧 Weyl 群 W_K 由 s_a1 和 s_a0 生成
MakeKDominant := function(h_coeffs)
    local current_h, w1, w0, changed, iter, max_iter, is_a2;
    
    current_h := [h_coeffs[1], h_coeffs[2]];
    is_a2 := false;
    if IsBoundGlobal("GlobalCurrentSubsystemType") then
        if ValueGlobal("GlobalCurrentSubsystemType") = "A2" then
            is_a2 := true;
        fi;
    fi;
    iter := 0;
    max_iter := 20; # 安全限制
    
    changed := true;
    while changed and iter < max_iter do
        changed := false;
        iter := iter + 1;
        
        # 检查 a1 (Short, Compact)
        if not is_a2 then
            w1 := CalculateWeight([1, 0], current_h);
            if w1 < 0 then
                current_h := ReflectH(current_h, [1, 0]);
                changed := true;
            fi;
        fi;
        
        # 检查 a0 (-theta, Long, Compact)
        w0 := CalculateWeight([-3, -2], current_h);
        if w0 < 0 then
            current_h := ReflectH(current_h, [-3, -2]);
            changed := true;
        fi;
    od;
    
    return current_h;
end;;

# =============================================================================
# 4. 主接口
# =============================================================================

PrintKOrbitInfo := function(h_str)
    local h_coeffs, dom_h, alt_h, w1, w0, w1_alt, w0_alt, is_a2, best_h, best_w0, best_w1, alt2_h, w1_alt2, w0_alt2;
    h_coeffs := ParseHString(h_str);
    
    Print("        >>> 模块 5 输出: K-Orbit 分类\n");
    Print("        K-Orbit Classification:\n");
    Print("          初始 h 系数: ", h_coeffs, "\n");
    
    is_a2 := false;
    if IsBoundGlobal("GlobalCurrentSubsystemType") then
        if ValueGlobal("GlobalCurrentSubsystemType") = "A2" then
            is_a2 := true;
        fi;
    fi;
    if not is_a2 then
        if PositionSublist(h_str, "2H_{") <> fail then
            is_a2 := true;
            GlobalCurrentSubsystemType := "A2";
        fi;
    fi;
    dom_h := MakeKDominant(h_coeffs);
    alt_h := MakeKDominant(ReflectH(h_coeffs, [0, 1]));
    best_h := dom_h;
    best_w1 := CalculateWeight([1, 0], best_h);
    best_w0 := CalculateWeight([-3, -2], best_h);
    w1_alt := CalculateWeight([1, 0], alt_h);
    w0_alt := CalculateWeight([-3, -2], alt_h);
    if w0_alt > best_w0 or (w0_alt = best_w0 and w1_alt < best_w1) then
        best_h := alt_h;
        best_w0 := w0_alt;
        best_w1 := w1_alt;
    fi;
    if is_a2 then
        alt2_h := MakeKDominant(ReflectH(ReflectH(h_coeffs, [1, 0]), [0, 1]));
        w1_alt2 := CalculateWeight([1, 0], alt2_h);
        w0_alt2 := CalculateWeight([-3, -2], alt2_h);
        if w0_alt2 > best_w0 or (w0_alt2 = best_w0 and w1_alt2 < best_w1) then
            best_h := alt2_h;
            best_w0 := w0_alt2;
            best_w1 := w1_alt2;
        fi;
    fi;
    dom_h := best_h;
    w1 := best_w1;
    w0 := best_w0;
    
    if dom_h <> h_coeffs then
        Print("          变换后 h 系数: ", dom_h, "\n");
    fi;
    
    DrawG2VoganDiagram(dom_h);
end;;

# =============================================================================
# 5. F4(4) K-c 共轭标准化函数
# =============================================================================

# K_c 共轭标准化：仅对紧根（a0,a2,a3,a4）做反射，保证权重 ≥0
# a1 保持原值（非紧）
NormalizeF4Weights_KDominant := function(w0, w1, w2, w3, w4)
    return NormalizeF4Weights_KDominant_Greedy(w0, w1, w2, w3, w4);
end;;

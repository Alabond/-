# G2 的简单根内积矩阵。
# 这里采用简单根基下的 Gram 矩阵，用于长度、角度与 Cartan 数据计算。
G2_InnerProdMatrix := [[2, -3], [-3, 6]];;

# 计算两个 G2 简单根坐标向量的内积。
# 输入是形如 [m1, m2] 的简单根系数向量，输出对应双线性型的值。
InnerProduct := function(v1, v2)
    return v1[1] * G2_InnerProdMatrix[1][1] * v2[1] +
           v1[1] * G2_InnerProdMatrix[1][2] * v2[2] +
           v1[2] * G2_InnerProdMatrix[2][1] * v2[1] +
           v1[2] * G2_InnerProdMatrix[2][2] * v2[2];
end;;

G2_SLAData := function()
    local L, R, S, Sinv, pos, roots;
    if IsBoundGlobal("GlobalG2SLAData") then
        return ValueGlobal("GlobalG2SLAData");
    fi;
    if not (IsBoundGlobal("LoadPackage") and LoadPackage("sla") = true) then
        Error("需要 sla 包");
    fi;
    L := SimpleLieAlgebra("G", 2, Rationals);
    R := RootSystem(L);
    S := SimpleRootsAsWeights(R);
    Sinv := S^-1;
    pos := List(PositiveRootsAsWeights(R), w -> List(w * Sinv, x -> Int(x)));
    roots := Concatenation(pos, List(pos, x -> -x));
    BindGlobal("GlobalG2SLAData", rec(all_roots := roots));
    return ValueGlobal("GlobalG2SLAData");
end;;

G2_AllRoots := function()
    return G2_SLAData().all_roots;
end;;

G2_WeylImages := function()
    local simple_roots, ReflectLocal, queue, imgs, seen, pair, w1, w2, i, nw1, nw2, key;
    if IsBoundGlobal("GlobalG2WeylImages") then
        return ValueGlobal("GlobalG2WeylImages");
    fi;
    simple_roots := [[1, 0], [0, 1]];
    ReflectLocal := function(v, alpha)
        local coeff;
        coeff := 2 * InnerProduct(v, alpha) / InnerProduct(alpha, alpha);
        return v - coeff * alpha;
    end;
    queue := [[[1, 0], [0, 1]]];
    imgs := [];
    seen := [];
    while Length(queue) > 0 do
        pair := queue[Length(queue)];
        Unbind(queue[Length(queue)]);
        w1 := pair[1];
        w2 := pair[2];
        key := Concatenation(String(w1), "|", String(w2));
        if Position(seen, key) = fail then
            Add(seen, key);
            Add(imgs, [w1, w2]);
            for i in [1..2] do
                nw1 := ReflectLocal(w1, simple_roots[i]);
                nw2 := ReflectLocal(w2, simple_roots[i]);
                Add(queue, [nw1, nw2]);
            od;
        fi;
    od;
    BindGlobal("GlobalG2WeylImages", imgs);
    return imgs;
end;;

# 判断一个根在当前 G2 染色规则下是否为黑点。
# 这里采用 a1 系数奇偶性作为非紧/紧的简化判据。
IsRootBlack := function(root)
    local n2;
    n2 := root[2];
    if n2 < 0 then
        n2 := -n2;
    fi;
    return (n2 mod 2) <> 0;
end;;

# 将简单根坐标向量格式化为 a1、a2 线性组合字符串。
# 例如 [2, -1] 会被格式化为 "2a1-a2"。
FormatRoot := function(v)
    local s;
    s := "";
    if v[1] <> 0 then
        if v[1] = 1 then
            s := Concatenation(s, "a1");
        elif v[1] = -1 then
            s := Concatenation(s, "-a1");
        else
            s := Concatenation(s, String(v[1]), "a1");
        fi;
    fi;
    if v[2] <> 0 then
        if v[1] <> 0 and v[2] > 0 then
            s := Concatenation(s, "+");
        fi;
        if v[2] = 1 then
            s := Concatenation(s, "a2");
        elif v[2] = -1 then
            s := Concatenation(s, "-a2");
        else
            s := Concatenation(s, String(v[2]), "a2");
        fi;
    fi;
    if s = "" then
        s := "0";
    fi;
    return s;
end;;

# 枚举、筛选并返回 G2 子系统数据。
# 该函数支持三种模式：
# 1) SCAN 模式：按黑白点数量和模板约束扫描候选子系统；
# 2) custom_roots 模式：直接对用户给定根做分类与着色；
# 3) 默认模式：沿 Weyl 轨道生成标准 G2 扩展子系统候选。
# 返回值中的每一项都包含根列表、颜色列表、连通分量、黑白点数量和类型标签。
FindG2Subsystems := function(custom_roots, custom_colors)
    local WeylImages, idx, w_alpha1, w_alpha2, w_theta, neg_w_theta,
          subsystem, results, black_nodes, white_nodes, root_strs, root_vecs, i, color, root_v,
          type_info, colors_list, c1, c2, target_black, target_white, all_roots, idx_count,
          r1, r2, ip, j, desired_type, desired_long_long, t1, t2, len1, len2, seen, key,
          l1, l2, tmp, tmpc, type_str, parser, ParseRootStringLocal, mid_color,
          pos_roots, pos_root_strs, pos_colors, pos_labels_local, i_local, j_local;

    results := [];

    if IsRecord(custom_roots) and IsBound(custom_roots.mode) and custom_roots.mode = "SCAN" then
        target_black := custom_roots.target_black;
        target_white := custom_roots.target_white;
        desired_type := fail;
        desired_long_long := false;
        if IsBound(custom_roots.template) and Length(custom_roots.template) = 2 then
            t1 := custom_roots.template[1];
            t2 := custom_roots.template[2];
            len1 := InnerProduct(t1, t1);
            len2 := InnerProduct(t2, t2);
            if target_black = 1 and target_white = 1 and len1 = 6 and len2 = 6 then
                desired_type := "A2";
                desired_long_long := true;
            elif target_black = 2 and target_white = 0 then
                desired_type := "A1 + A1";
            fi;
        fi;

        if desired_type = "A2" and desired_long_long then
            WeylImages := G2_WeylImages();
            idx_count := 0;
            seen := [];
            for idx in [1..Length(WeylImages)] do
                w_alpha1 := WeylImages[idx][1];
                w_alpha2 := WeylImages[idx][2];
                w_theta := 3 * w_alpha1 + 2 * w_alpha2;
                neg_w_theta := -w_theta;
                c1 := 0;
                if IsRootBlack(w_alpha2) then
                    c1 := 1;
                fi;
                c2 := 0;
                if IsRootBlack(neg_w_theta) then
                    c2 := 1;
                fi;
                if (c1 + c2) <> target_black or ((1 - c1) + (1 - c2)) <> target_white then
                    continue;
                fi;
                type_info := ClassifyRootSystem([w_alpha2, neg_w_theta], InnerProduct);
                if type_info.type_str <> desired_type then
                    continue;
                fi;
                r1 := neg_w_theta;
                r2 := w_alpha2;
                root_strs := [FormatRoot(r1), FormatRoot(r2)];
                key := Concatenation(root_strs[1], "|", root_strs[2]);
                if Position(seen, key) <> fail then
                    continue;
                fi;
                Add(seen, key);
                idx_count := idx_count + 1;
                colors_list := [c2, c1];
                mid_color := 0;
                if IsRootBlack(w_alpha1) then
                    mid_color := 1;
                fi;
                Add(results, rec(
                    idx := idx_count,
                    roots_list := root_strs,
                    all_roots_list := [FormatRoot(neg_w_theta), FormatRoot(w_alpha1), FormatRoot(w_alpha2)],
                    all_colors_list := [c2, mid_color, c1],
                    pos_labels := ["a0", "a2"],
                    all_pos_labels := ["a0", "a1", "a2"],
                    colors_list := colors_list,
                    components := type_info.components,
                    black_nodes := c1 + c2,
                    white_nodes := (1 - c1) + (1 - c2),
                    type := type_info.type_str,
                    w_alpha2 := FormatRoot(w_alpha2),
                    neg_w_theta := FormatRoot(neg_w_theta)
                ));
            od;
            return results;
        fi;
        
        if desired_type = "A1 + A1" then
            WeylImages := G2_WeylImages();
            idx_count := 0;
            seen := [];
            for idx in [1..Length(WeylImages)] do
                w_alpha1 := WeylImages[idx][1];
                w_alpha2 := WeylImages[idx][2];
                w_theta := 3 * w_alpha1 + 2 * w_alpha2;
                neg_w_theta := -w_theta;
                pos_roots := [neg_w_theta, w_alpha1, w_alpha2];
                pos_root_strs := List(pos_roots, r -> FormatRoot(r));
                pos_labels_local := ["a0", "a1", "a2"];
                pos_colors := [];
                for i_local in [1..3] do
                    if IsRootBlack(pos_roots[i_local]) then
                        Add(pos_colors, 1);
                    else
                        Add(pos_colors, 0);
                    fi;
                od;
                for i_local in [1..2] do
                    for j_local in [i_local + 1..3] do
                        if pos_colors[i_local] <> 1 or pos_colors[j_local] <> 1 then
                            continue;
                        fi;
                        if InnerProduct(pos_roots[i_local], pos_roots[j_local]) <> 0 then
                            continue;
                        fi;
                        type_info := ClassifyRootSystem([pos_roots[i_local], pos_roots[j_local]], InnerProduct);
                        if type_info.type_str <> desired_type then
                            continue;
                        fi;
                        root_strs := [pos_root_strs[i_local], pos_root_strs[j_local]];
                        key := Concatenation(root_strs[1], "|", root_strs[2]);
                        if Position(seen, key) <> fail then
                            continue;
                        fi;
                        Add(seen, key);
                        idx_count := idx_count + 1;
                        Add(results, rec(
                            idx := idx_count,
                            roots_list := root_strs,
                            all_roots_list := pos_root_strs,
                            all_colors_list := ShallowCopy(pos_colors),
                            pos_labels := [pos_labels_local[i_local], pos_labels_local[j_local]],
                            all_pos_labels := pos_labels_local,
                            colors_list := [1, 1],
                            components := type_info.components,
                            black_nodes := 2,
                            white_nodes := 0,
                            type := type_info.type_str,
                            w_alpha2 := root_strs[1],
                            neg_w_theta := root_strs[2]
                        ));
                    od;
                od;
            od;
            return results;
        fi;

        all_roots := G2_AllRoots();
        idx_count := 0;
        for i in [1..Length(all_roots)] do
            for j in [i + 1..Length(all_roots)] do
                r1 := all_roots[i];
                r2 := all_roots[j];
                c1 := 0;
                if IsRootBlack(r1) then
                    c1 := 1;
                fi;
                c2 := 0;
                if IsRootBlack(r2) then
                    c2 := 1;
                fi;
                if (c1 + c2) <> target_black or ((1 - c1) + (1 - c2)) <> target_white then
                    continue;
                fi;
                type_info := ClassifyRootSystem([r1, r2], InnerProduct);
                if desired_type <> fail then
                    if type_info.type_str <> desired_type then
                        continue;
                    fi;
                    if desired_long_long then
                        if InnerProduct(r1, r1) <> 6 or InnerProduct(r2, r2) <> 6 then
                            continue;
                        fi;
                    fi;
                fi;
                if desired_type = "A1 + A1" then
                    ip := InnerProduct(r1, r2);
                    if ip <> 0 then
                        continue;
                    fi;
                fi;
                l1 := InnerProduct(r1, r1);
                l2 := InnerProduct(r2, r2);
                if l1 < l2 then
                    tmp := r1;
                    r1 := r2;
                    r2 := tmp;
                    tmpc := c1;
                    c1 := c2;
                    c2 := tmpc;
                fi;
                idx_count := idx_count + 1;
                root_strs := [FormatRoot(r1), FormatRoot(r2)];
                colors_list := [c1, c2];
                Add(results, rec(
                    idx := idx_count,
                    roots_list := root_strs,
                    colors_list := colors_list,
                    components := type_info.components,
                    black_nodes := c1 + c2,
                    white_nodes := (1 - c1) + (1 - c2),
                    type := type_info.type_str,
                    w_alpha2 := root_strs[1],
                    neg_w_theta := root_strs[2]
                ));
            od;
        od;
        return results;
    fi;

    if custom_roots <> fail then
        ParseRootStringLocal := function(root_str)
            local v, coeff, idx_local, sign, num_str;
            v := [0, 0];
            root_str := Filtered(root_str, c -> c <> ' ');
            if root_str = "0" then
                return v;
            fi;
            idx_local := 1;
            while idx_local <= Length(root_str) do
                sign := 1;
                if root_str[idx_local] = '-' then
                    sign := -1;
                    idx_local := idx_local + 1;
                elif root_str[idx_local] = '+' then
                    idx_local := idx_local + 1;
                fi;
                num_str := "";
                while idx_local <= Length(root_str) and root_str[idx_local] in "0123456789" do
                    num_str := Concatenation(num_str, [root_str[idx_local]]);
                    idx_local := idx_local + 1;
                od;
                if num_str = "" then
                    coeff := 1;
                else
                    coeff := Int(num_str);
                fi;
                coeff := coeff * sign;
                if idx_local <= Length(root_str) and root_str[idx_local] = 'a' then
                    idx_local := idx_local + 1;
                    if idx_local <= Length(root_str) and root_str[idx_local] in "12" then
                        if root_str[idx_local] = '1' then
                            v[1] := v[1] + coeff;
                        else
                            v[2] := v[2] + coeff;
                        fi;
                        idx_local := idx_local + 1;
                    fi;
                else
                    break;
                fi;
            od;
            return v;
        end;
        parser := ParseRootStringLocal;
        if IsBoundGlobal("ParseRootString") then
            parser := ValueGlobal("ParseRootString");
        fi;
        root_vecs := [];
        for i in [1..Length(custom_roots)] do
            Add(root_vecs, parser(custom_roots[i]));
        od;
        black_nodes := 0;
        white_nodes := 0;
        if custom_colors <> fail then
            if Length(custom_colors) <> Length(custom_roots) then
                Error("颜色数量与根数量不匹配");
            fi;
            for color in custom_colors do
                if color = 1 then
                    black_nodes := black_nodes + 1;
                else
                    white_nodes := white_nodes + 1;
                fi;
            od;
        else
            for root_v in root_vecs do
                if IsRootBlack(root_v) then
                    black_nodes := black_nodes + 1;
                else
                    white_nodes := white_nodes + 1;
                fi;
            od;
        fi;
        type_info := ClassifyRootSystem(root_vecs, InnerProduct);
        type_str := type_info.type_str;
        if custom_colors = fail then
            colors_list := [];
            for root_v in root_vecs do
                if IsRootBlack(root_v) then
                    Add(colors_list, 1);
                else
                    Add(colors_list, 0);
                fi;
            od;
        else
            colors_list := custom_colors;
        fi;
        Add(results, rec(
            idx := 0,
            w_alpha2 := custom_roots[1],
            neg_w_theta := custom_roots[Length(custom_roots)],
            roots_list := custom_roots,
            colors_list := colors_list,
            components := type_info.components,
            black_nodes := black_nodes,
            white_nodes := white_nodes,
            type := type_str
        ));
        return results;
    fi;

    WeylImages := G2_WeylImages();

    for idx in [1..Length(WeylImages)] do
        w_alpha1 := WeylImages[idx][1];
        w_alpha2 := WeylImages[idx][2];
        w_theta := 3 * w_alpha1 + 2 * w_alpha2;
        neg_w_theta := -w_theta;
        black_nodes := 0;
        white_nodes := 0;
        c1 := 0;
        if IsRootBlack(w_alpha2) then
            black_nodes := black_nodes + 1;
            c1 := 1;
        else
            white_nodes := white_nodes + 1;
        fi;
        c2 := 0;
        if IsRootBlack(neg_w_theta) then
            black_nodes := black_nodes + 1;
            c2 := 1;
        else
            white_nodes := white_nodes + 1;
        fi;
        if c1 <> c2 then
            subsystem := [w_alpha2, neg_w_theta];
            type_info := ClassifyRootSystem(subsystem, InnerProduct);
            mid_color := 0;
            if IsRootBlack(w_alpha1) then
                mid_color := 1;
            fi;
            Add(results, rec(
                idx := idx,
                w_alpha2 := FormatRoot(w_alpha2),
                neg_w_theta := FormatRoot(neg_w_theta),
                roots_list := [FormatRoot(w_alpha2), FormatRoot(neg_w_theta)],
                all_roots_list := [FormatRoot(neg_w_theta), FormatRoot(w_alpha1), FormatRoot(w_alpha2)],
                all_colors_list := [c2, mid_color, c1],
                pos_labels := ["a2", "a0"],
                all_pos_labels := ["a0", "a1", "a2"],
                colors_list := [c1, c2],
                components := type_info.components,
                black_nodes := black_nodes,
                white_nodes := white_nodes,
                type := type_info.type_str
            ));
        fi;
    od;

    return results;
end;;

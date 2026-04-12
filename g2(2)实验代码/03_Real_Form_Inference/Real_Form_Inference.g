#############################################################################
# 文件: Real_Form_Inference.g
# 描述: 根据 Dynkin 图类型和 Vogan 图颜色推断实李代数类型
# 说明: 模块2规则入口文件。
# 函数用途:
#   - InferRealForm(type, rank, black_nodes, white_nodes):
#       按类型与黑白计数输出实形式记录 rec(type, params, desc)。
# 备注:
#   - 当前映射为项目实验导向的启发式规则，强调“流程可比对与可复现”。
#############################################################################

# 推断实李代数类型
# type: "A", "B", "C", "D"
# rank: int
# black_nodes: int (非紧根数量)
# white_nodes: int (紧根数量)
# total_nodes = black_nodes + white_nodes = rank
#
# 返回: rec(real_form_type, params)
# 
# 简化规则（面向本项目的 G2 子系统场景）：
# - A: 全黑→sl(n+1,ℝ)；全白→紧（su(n+1)）；混合→su(p,q) 的启发式映射
# - C: 全黑→sp(n,ℝ)；否则→sp(p,q)（占位参数）
# - B/D: 全黑→近似 split 形 so(n+1,n)；否则→so(p,q)（占位参数）

# 根据根系类型、秩与黑白点数量推断对应实形式。
# 当前规则是为 G2 子系统实验定制的启发式映射，输出记录包含类型编码、参数与文字说明。
InferRealForm := function(type, rank, black_nodes, white_nodes)
    local p, q;
    
    if type = "A" then
        if rank = 1 then
            if black_nodes = 1 then
                return rec(type := "sl_R", params := [1], desc := "sl(2, R)");
            fi;
            return rec(type := "compact", params := [2], desc := "su(2)");
        fi;
        if white_nodes = rank then
            return rec(type := "compact", params := [rank+1], desc := Concatenation("su(", String(rank+1), ")"));
        elif black_nodes = rank then
            return rec(type := "sl_R", params := [rank], desc := Concatenation("sl(", String(rank+1), ", R)"));
        elif rank = 3 and ((black_nodes = 1 and white_nodes = 2) or (black_nodes = 2 and white_nodes = 1)) then
            return rec(type := "su_pq", params := [2, 2], desc := "su(2,2)");
        else
            p := white_nodes + 1;
            q := black_nodes;
            if q = 0 then q := 1; p := rank; fi; 
            return rec(type := "su_pq", params := [p, q], desc := Concatenation("su(", String(p), ",", String(q), ")"));
        fi;
        
    elif type = "C" then
        if black_nodes = rank then
            return rec(type := "sp_R", params := [rank], desc := "sp(n, R)");
        else
            return rec(type := "sp_pq", params := [1, 1], desc := "sp(p, q)");
        fi;
        
    elif type = "B" or type = "D" then
        if black_nodes = rank then
            # Split form for B/D is so(n, n+1) or so(n, n)
            # so(p, q) split: p ~ q
            return rec(type := "so_pq", params := [rank+1, rank], desc := "so(split)");
        else
            return rec(type := "so_pq", params := [rank+1, 1], desc := "so(p, q)");
        fi;
    fi;
    
    return fail;
end;;

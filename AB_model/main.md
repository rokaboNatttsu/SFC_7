- [1. 問題とTODO](#1-問題とtodo)
- [2. 推奨実行環境](#2-推奨実行環境)
- [3. エージェント番号](#3-エージェント番号)
- [4. モデルの式とアルゴリズム](#4-モデルの式とアルゴリズム)
- [5. 拡張と改変の方向性](#5-拡張と改変の方向性)

# 1. 問題とTODO

優先順に番号を振っている

1. jlファイルの変更をmain.mdに反映させる・説明追加

# 2. 推奨実行環境
- RAM：最低8GB以上、推奨16GB以上
  - ほとんどミニマムな条件
  - J, O, N, TIME,  = 1000, 10, 3, 200
  - Oin = 1
  - OBuffer = O + TIME*Oin
  - でも、実行のために4GB程度使っていた

# 3. エージェント番号
|      | インデックス | 総数 |
| :--- | :----------: | :--: |
| 家計 |     $j$      | $J$  |
| 企業 |     $o$      | $O$  |
| 銀行 |     $n$      | $N$  |

変数の右上につけて、どのエージェントとどのエージェントの取引もしくは貸借関係なのかを区別する

# 4. モデルの式とアルゴリズム

変数の定義を兼ねた、会計的整合性を表す表は、matrix.mdに記載している。

randn() : 標準化された正規分布で確率分布が表される乱数からのサンプル

計算順に式を一覧する

1. $p^o=\nu_2(1+\nu_1+\nu_3\frac{i^o_{-1}}{u^T \gamma_1 k^o_{-1}})\frac{\sum_j W^{j,o}_{-1}+T^o_{v-1}+T^o_{c-1}+r_L\sum_n L^{o,n}_{f-1}+\delta k^o_{-1}}{u^T \gamma_1 k^o_{-1}} + (1-\nu_2)\bar{p}$
   1. ここで、 $\bar{p}=\frac{\sum_o p^o_{-1}(\sum_j c^{o,j}_{-1} + g^o_{-1})}{\sum_o(\sum_j c^{o,j}_{-1}+g^o_{-1})}$
2. 就業・失業判定、賃金率と賃金の決定のアルゴリズム
   1. 前期の労働稼働率が１より大きいときは、労働稼働率が１になる水準の求人を出す
   2. 前期就業していた人の一定割合はその企業での就業を続行
   3. 前期就業していた人のうち、勤め先が倒産した人は、一定確立で即座に求人に応募する
   4. 前期就業していて勤め先が倒産していない人のうち、一定割合がその仕事をやめ、一定確立で即座に求人に応募
   5. 前期失業していた人の中からランダムで、起業する
   6. 前期失業していて企業もしない人は、求人に応募する。応募先は求人数に比例してランダムに選ぶ
   7. 求人数が満たされるまでランダムで応募してきた人を雇う
   8. 前期就業していた人に賃金を支払い。賃金の支払いは前期就業していた企業 $o'$ から行われる。金額は、その企業の前期の労働稼働率 $v^{o'}_{-1}$ と家計の要求賃金率 $w^j_{-1}$ の積 $W^{j,o'}=v^{o'}_{-1}w^j_{-1}$
   9. 従業員がいなくなった企業を倒産させる
   10. 今期就業中の家計の要求賃金率 $w^j$ を前期の値 $w^j_{-1}$ より上げる。 $w^j=w^j_{-1}(1+\zeta_2 \times| randn() |)$
   11. 今期失業中の家計の要求賃金率 $w^j$ を前期の値 $w^j_{-1}$ より下げる。 $w^j=w^j_{-1}(1-\zeta_3 \times | randn() |)$
3.  $T_i^j=\tau_1 \sum_j (W_{-1}^{j,o}+IT^j_{h-1}+P_{h-1}^{j,o}+S_{-1}^{j,n})$
4.  $T_a^j=\tau_2(\sum_o E_{h-1}^{o,j}+\sum_n F^{n,j}+\sum_n M_{h-1}^{n,j}-\sum_n L_{h-1}^{j,n})$
5. $T_v^o=\tau_3(\sum_j C_{-1}^{o,j}+I_{-1}^o+G_{-1}^o)$
6. $T_c^o=\max(0, \tau_4 P^o_{-1})$
7. 名目政府支出の需要のアルゴリズム
   1. $G_{sum}=G_0(1+\beta_1)^{(t-1)}$
   2. $x^o=\max[1, x_{-1}^o \{1-\beta_2+\beta_3 randn()\}]$
   3. $z^o=\max(0, x^o-\beta_4) (k^o_{t-1})^{\beta_5}$
   4. $G^{o,D}=\frac{z^o G_{sum}}{\sum z^o}$
   5. 下限付きべき分布みたいなものを作って(xがこれにあたる)、一定の値を差し引いて保有する実質資本の$\beta_5$乗した(zがこれに当たる)後に、政府支出総額に合わせて全体的に値を大きくする。
8.  $g^{o,D}=\frac{G^{o,D}}{p^o}$
9.  家計の消費需要
    1. $Cs^{o}=\max \{ \alpha_4\frac{\sum_j\sum_o C^{o,j}_{-1}}{J}, \alpha_1(\sum_o W_{-1}^{o,j}-T_a^j-T_i^j-r_L \sum_n L_{h-1}^{j,n}+IT^j_{h-1}+\sum_o P_{h-1}^{o,j}+\sum_n S_{-1}^{n,j}) + \alpha_2(\sum_o E_{h-1}^{o,j}+\sum_n M_{h-1}^{n,j}-\sum_n L_{h-1}^{n,j}) \}$
10. 消費財の購入先選択アルゴリズム
    1.  前回消費財を買った企業$o$から別の企業に乗り換える確率は $rand() < \alpha_3 \frac{p^{o}-\overline{p}}{\overline{p}}$ とする。 $\overline{p}$ は平均価格で、前期の家計の消費と政府支出の名目値の合計を、前期の家計の消費と政府支出の実質値で割った値。消費の乗り換え先は「前期の実質生産量($\sum_j c_{-1}^{o,j}+g_{-1}^o$)」に比例する確率でランダムに決める
    2.  今期の消費財購入先の企業を$o'$と書くことにすると、今期の消費需要は、 $c^{o',j,D}=\frac{Cs^{o'}}{p^{o'}}$ で決める。 $o'$ 以外のすべての企業の商品は消費しない( $c^{o,j}=0 \ \ \ \{ o \mid o \neq o'\}$ )
11. $i^{o,D}=\max[0, \min \{ \delta k_{-1}^o + (u^o-u^T)\gamma_1 k_{-1}^o, \gamma_2\frac{\sum_n M_{f-1}^{n,o}-\sum_n L_{f-1}^{o,n}}{p^o}\}]$
12. $u^{o,D} = \frac{\sum_j c^{o,j,D}+i^{o,D}+g^{o,D}}{\gamma_1 k_{-1}^o}$
13. 資本稼働率 $u^o$ が１以下になるように、需要が多すぎるときは供給量を減らす。
    1.  $u^o = \frac{u^{o,D}}{\max(1, u^{o,D})}$
    2.  $G^o = \frac{G^{o,D}}{\max(1, G^{o,D})}$
    3.  $g^o = \frac{g^{o,D}}{\max(1, g^{o,D})}$
    4.  $c^{o,j} = \frac{c^{o,j,D}}{\max(1, c^{o,j,D})}$
    5.  $i^o = \frac{i^{o,D}}{\max(1, i^{o,D})}$
14. $C^{o,j}=p^o c^{o,j}$
15. $I^o=p^o i^o$
16. $A^o=A_{-1}^o(1 + \mu_1 + \mu_2 \frac{i_{-1}^o}{k_{-1}^o})$
    1.  技術水準上昇率が投資の一次関数として決まる。技術水準は、労働者一人時間あたりが扱う資本の量で表される。技術水準が上がるほど資本の量が増えるのでは？という発想。「情報と秩序」というタイトルの本からインスピレーションを受けている
17. $k^o=(1-\delta)k_{-1}^o+i^o$
18. $K^o=p^o k^o$
19. 政府から家計への所得移転
    1.  失業しているとき $IT^j_h=(1+\chi_1+\chi_2) \frac{\sum_j\sum_o C^{o,j}_{-1}}{J}$
    2.  就業しているとき $IT^j_h=(1+\chi_1) \frac{\sum_j\sum_o C^{o,j}_{-1}}{J}$
21. $P^o=\sum_j C^{o,j}+G^o+I^o-\sum_j W^{j,o}-T_c^o-T_v^o-r_L \sum_n L_{f-1}^{o,n}$
22. $P_h^{j,o}=\max\{0, \theta_1(P^o-I^o)+\theta_2(\sum_n M_{f-1}^{o,n}-\sum_n L_{f-1}^{o,n})\}\frac{e_{h-1}^{o,j}}{e_{-1}^o}$
23. $P_b^{n,o}=\max\{0, \theta_1(P^o-I^o)+\theta_2(\sum_n M_{f-1}^{o,n}-\sum_n L_{f-1}^{o,n})\}\frac{e_{b-1}^{o,n}}{e_{-1}^o}$
24. $P_f^o=P^o-\sum_j P_h^{j,o}-\sum_n P_b^{n,o}$
25. $S^{j,n}=\max[\{\theta_1(r_L L_{-1}^n+\sum_o P_b^{n,o})+\theta_2 \sum_o E_{b-1}^{o,n}\}\frac{f_{h-1}^{j,n}}{f_{-1}^n}]$
- 
- ここから説明更新
- 
- $NL_h^j=-\sum_o C^{o,j}+\sum_o W^{j,o}-T_i^j-T_a^j-r_L \sum_n L_{h-1}^{j,n}+\sum_o P_h^{j,o}+\sum_n S^{j,n}$
- $NL_f^o=-I^o+P_f^o$
- $NL_b^n=r_L L_{-1}^n+\sum_o P_b^{o,n}-\sum_j S^{j,n}$
- $NL_g=-\sum_o G^o+\sum_j T_i^j+\sum_j T_a^j+\sum_o T_v^o+\sum_o T_c^o$
- $\sum_o L_h^{j,n}=\epsilon_1 NL_h^j+\epsilon_2 \sum_o C^{o,j}$
  - 家計は預金口座を持つ銀行を借入先に選ぶ
    - $L_h^{j,n}=$
    - 本当は銀行ごとに信用スコアを計算して金利が安いところから借り入れようとするとか、既存の借入先を維持する傾向を持つとか、借入先を2つ以上にすることもあるとか、そういう効果を入れたいが、モデルの複雑さを抑えるためのアドホックな仮定として導入することにする。
- $\Delta L_h^{j,n}=L_h^{j,n}-L_{h-1}^{j,n}$
- $\sum_n \Delta L_f^{o,n}=\max\{-\sum_n L_{f-1}^{o,n}, (\lambda_3 + \lambda_4(\frac{P^o - P_f^o}{\sum_j E_{h-1}^{o,j} + \sum_n E_{b-1}^{o,n}} - r_L))(I^o+\sum_j W^{j,o}+T_v^o+T_c^o+r_L \sum_n L_{f-1}^{j,n} - \phi \sum_n M_{f-1}^{n,j})\}$
  - もっとリアルな貸付金水準の行動方程式はないか？要調査
  - $L_h^{j,n}$ と同じく、企業が預金口座を持つ銀行から借り入れる
- $L_f^{o,n}=L_{f-1}^{o,n}+\Delta L_f^{o,n}$
- $L^n=\sum_j L_h^{j,n}+\sum_o L_f^{o,n}$
- $\Delta L^n=L^n-L_{-1}^n$
- $\Delta e^o=\frac{1}{p_{-1}^o}\{1-\lambda_3 - \lambda_4(\frac{P^o - P_f^o}{\sum_j E_{h-1}^{o,j} + \sum_n E_{b-1}^{o,n}} - r_L)\}(I^o+\sum_j W^{j,o}+T_v^o+T_c^o+r_L \sum_n L_{f-1}^{o,n} - \phi \sum_n M_{f-1}^{n,o})$
- $E^o=E_{-1}^o+p_{e-1}^o\Delta e^o$
- $e^o=e_{-1}^o+\Delta e^o$
- 家計はどうやって株式保有額を決めるのか
  - どこの株を保有するかの意思決定
    - 利益、純金融資産、配当率　で重み付けされたランダムサンプリングと、前期保有していた銘柄を持ち続ける傾向の合成
  - $E_h^j=$
  - $F_h^j=$
- 銀行の企業株式保有割合を決める
  - 利益、純金融資産、配当率　で重み付けされた保有割合に、摂動を加え、期首の保有割合の影響も残す
  - $E_b^{o,n}=$
- $p_e^o=\frac{\sum_j E_h^{o,j}+\sum_n E_b^{o,n}}{e^o}$
- $e_h^{o,j}=\frac{E_h^{o,j}}{p_e^o}$
- $e_b^{o,n}=\frac{E_b^{o,n}}{p_e^o}$
- $\Delta e_h^{o,j}=e_h^{o,j}-e_{h-1}^{o,j}$
- $\Delta e_b^{o,n}=e_h^{o,n}-e_{b-1}^{o,n}$
- 家計の預金の決定
  - $\sum_n \Delta M_h^{n,j}=NL_h^j-\sum_o(p_{e-1}^o\Delta e_h^{o,j})+\sum_n \Delta L_h^{j,n}$
  - 一定確率で、預金の保有先を変更する。変更先は、純資産に比例する確率で選ぶ
  - $M_h=M_{h-1}+\Delta M_h$
- 企業の預金の決定
  - $\sum_n \Delta M_f^{n,o}=NL_f^o+\sum_n \Delta L_f^{o,n}+p_{e-1}^o\Delta e^o$
  - 一定確率で、預金の保有先を変更する。変更先は、純資産に比例する確率で選ぶ
  - $M_f^{n,o}=M_{f-1}^{n,o}+\Delta M_f^{n,o}$
- $M^n=\sum_j M_h^{n,j}+\sum_o M_f^{n,o}$
- $\Delta M^n=M^n-M_{-1}^n$
- $\Delta H^n = NL_b^n - \sum_o p_{e-1}^o\Delta e_b^{o,n} + \Delta M^n - \Delta L^n$
- $H^n = H_{-1}^n + \Delta H^n$
- $NW_h^j=\sum_n M_h^{n,j}-\sum_n L_h^{j,n}+\sum_o E_h^{o,j}$
- $NW_f^o=K^o+\sum_n M_f^{n,o}-\sum_n L_f^{o,n}-E^o$
- $NW_b^n=-M^n+L^n+\sum_o E_b^{o,n}+H^n$
- $NW_g=-\sum_n H^n$
- $D_E=-\sum_o E^o+\sum_o\sum_j E_h^{o,j}+\sum_o\sum_j E_b^{o,j}$
- $D_F=\sum_j\sum_n F_h^{j,n}-\sum_n F^n$
- $NW=\sum_j NW_h^j+\sum_o NW_f+\sum_n NW_b+NW_g$
- $v^o=\frac{u^o k^o}{A^o \sum_j (W^{j,o}>0)}$
- 
- 倒産と新規参入のアルゴリズムも追加する必要がある
- 倒産
  - 純金融資産が0以下になると倒産する
  - 

隠された恒等式として

- 

を使ってモデルの会計的一貫性を確認する


# 5. 拡張と改変の方向性
- $r_{Lf}=(1-\lambda_{rLf})r_{Lf-1}+\lambda_{rLf}(1+\epsilon_3+\frac{\Delta Z^n}{\sum_o L_{f-1}^{o,n}})$ で企業向けの適応的な金利を導入することを検討
- 起業する企業の資本が何もないところから出てくるのではなく、既存企業の生産した資本を買う形式にする
- 政府や銀行が家計を雇う
- 国債の導入
- 銀行への資本規制の追加
- 金利の内生化
- 預金/借入金の多対多化
- 貸し渋り、貸し剥がしのような与信行動を追加する
- 失業者への社会保障をモデル含める
- など
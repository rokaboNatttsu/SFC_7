- [1. 集計的モデル](#1-集計的モデル)
  - [1.1. 会計的整合性を表す恒等式](#11-会計的整合性を表す恒等式)
    - [1.1.1. 取引フロー表](#111-取引フロー表)
    - [1.1.2. バランスシート表](#112-バランスシート表)
    - [1.1.3. 負債及び資産の増減と、純資産の増減の、整合性を表す表(Full-integration matrix)](#113-負債及び資産の増減と純資産の増減の整合性を表す表full-integration-matrix)
    - [1.1.4. ストックとフローの関係の整合性を示す表](#114-ストックとフローの関係の整合性を示す表)
    - [1.1.5. 表で示されない恒等式](#115-表で示されない恒等式)
  - [1.2. モデルの式一覧(計算する順)](#12-モデルの式一覧計算する順)
- [2. エージェント・ベース・モデル](#2-エージェントベースモデル)
  - [2.1. エージェント番号](#21-エージェント番号)
  - [2.2. 部門集計レベルの会計恒等式群](#22-部門集計レベルの会計恒等式群)
    - [2.2.1. 取引フロー表](#221-取引フロー表)
    - [2.2.2. バランスシート表](#222-バランスシート表)
    - [2.2.3. 負債及び資産の増減と、純資産の増減の、整合性を表す表(Full-integration matrix)](#223-負債及び資産の増減と純資産の増減の整合性を表す表full-integration-matrix)
  - [2.3. エージェントレベルの会計恒等式群](#23-エージェントレベルの会計恒等式群)
    - [2.3.1. フローの一貫性(貸借対照表及びキャッシュ・フロー計算書のようなもの)の表](#231-フローの一貫性貸借対照表及びキャッシュフロー計算書のようなものの表)
    - [2.3.2. ストックの一貫性(バランスシート)の表](#232-ストックの一貫性バランスシートの表)
    - [2.3.3. 負債及び資産の増減と、純資産の増減の、整合性を表す表(Full-integration matrix)](#233-負債及び資産の増減と純資産の増減の整合性を表す表full-integration-matrix)
    - [2.3.4. ストックとフローの関係の整合性を示す表](#234-ストックとフローの関係の整合性を示す表)
  - [2.4. 表で示されない恒等式](#24-表で示されない恒等式)
  - [2.5. モデルの式とアルゴリズム](#25-モデルの式とアルゴリズム)

# 1. 集計的モデル

いずれ作るAB-SFCモデルの基礎とする予定

simulation.ipynbで実装した

拡張する方向として検討すること

- 公債と公債利回りの追加
  - 公債利回りと民間への貸出金利の適応的連動 $r_L=(1-\lambda_r)r_{L-1}+\lambda_r \epsilon r_{GB}$ 
- 家計と投資家の分割
- 公務員の追加
- 政府の資本の追加
- 資本財生産企業と消費財生産企業の分割
- 海外の導入

## 1.1. 会計的整合性を表す恒等式
### 1.1.1. 取引フロー表
|                |         家計          |   企業(経常)   |     企業(資本)      |         銀行          |  統合政府   |  合計   |
| :------------- | :-------------------: | :------------: | :-----------------: | :-------------------: | :---------: | :-----: |
| 消費           |         $-C$          |      $+C$      |                     |                       |             |   $0$   |
| 政府支出       |                       |      $+G$      |                     |                       |    $-G$     |   $0$   |
| 投資           |                       |      $+I$      |        $-I$         |                       |             |   $0$   |
| 賃金           |         $+W$          |      $-W$      |                     |                       |             |   $0$   |
| 所得税         |        $-T_i$         |                |                     |                       |   $+T_i$    |   $0$   |
| 資産税         |        $-T_a$         |                |                     |                       |   $+T_a$    |   $0$   |
| 付加価値税     |                       |     $-T_v$     |                     |                       |   $+T_v$    |   $0$   |
| 法人税         |                       |     $-T_c$     |                     |                       |   $+T_c$    |   $0$   |
| 借入金利払     |    $-r_L L_{h-1}$     | $-r_L L_{f-1}$ |                     |     $+r_L L_{-1}$     |             |   $0$   |
| 企業利潤       |        $+P_h$         |      $-P$      |       $+P_f$        |        $+P_b$         |             |   $0$   |
| 銀行配当       |         $+S$          |                |                     |         $-S$          |             |   $0$   |
| [ 純貸出 ]     |      [ $NL_h$ ]       |                |     [ $NL_f$ ]      |      [ $NL_b$ ]       | [ $NL_g$ ]  | [ $0$ ] |
| 企業株式純発行 | $-p_{e-1} \Delta e_h$ |                | $+p_{e-1} \Delta e$ | $-p_{e-1} \Delta e_b$ |             |   $0$   |
| 預金の移動     |     $-\Delta M_h$     |                |    $-\Delta M_f$    |      $+\Delta M$      |             |   $0$   |
| 貸付金の移動   |     $+\Delta L_h$     |                |    $+\Delta L_f$    |      $-\Delta L$      |             |   $0$   |
| 現金の移動     |                       |                |                     |      $-\Delta H$      | $+\Delta H$ |   $0$   |
| 合計           |          $0$          |      $0$       |         $0$         |          $0$          |     $0$     |         |

### 1.1.2. バランスシート表
|          |  家計   |  企業   |  銀行   | 統合政府 | 合計  |
| :------- | :-----: | :-----: | :-----: | :------: | :---: |
| 資本     |         |  $+K$   |         |          | $+K$  |
| 預金     | $+M_h$  | $+M_f$  |  $-M$   |          |  $0$  |
| 貸付金   | $-L_h$  | $-L_f$  |  $+L$   |          |  $0$  |
| 企業株式 | $+E_h$  |  $-E$   | $+E_b$  |          | $+D$  |
| 現金     |         |         |  $+H$   |   $-H$   |  $0$  |
| 純資産   | $-NW_h$ | $-NW_f$ | $-NW_b$ | $-NW_g$  | $-NW$ |
| 合計     |   $0$   |   $0$   |   $0$   |   $0$    |  $0$  |

### 1.1.3. 負債及び資産の増減と、純資産の増減の、整合性を表す表(Full-integration matrix)
|                          |         家計          |        企業        |         銀行          |  統合政府   |         合計         |
| :----------------------- | :-------------------: | :----------------: | :-------------------: | :---------: | :------------------: |
| 期首純資産               |      $NW_{h-1}$       |     $NW_{f-1}$     |      $NW_{b-1}$       | $NW_{g-1}$  |      $NW_{-1}$       |
| 資本のキャピタルゲイン   |                       | $+\Delta p k_{-1}$ |                       |             |  $+\Delta p k_{-1}$  |
| 資本増減                 |                       |   $-p \Delta k$    |                       |             |    $-p \Delta k$     |
| 預金の増減               |     $+\Delta M_h$     |   $+\Delta M_f$    |      $-\Delta M$      |             |         $0$          |
| 貸付金の増減             |     $-\Delta L_h$     |   $-\Delta L_h$    |      $+\Delta L$      |             |         $0$          |
| 企業株式キャピタルゲイン | $+\Delta p_e e_{h-1}$ |                    | $+\Delta p_e e_{b-1}$ |             | $+\Delta p_e e_{-1}$ |
| 企業株式純発行           |   $+p_e \Delta e_h$   |  $+p_e \Delta e$   |   $+p_e \Delta e_b$   |             |         $0$          |
| 現金の増減               |                       |                    |      $+\Delta H$      | $-\Delta H$ |         $0$          |
| 純資産                   |        $-NW_h$        |      $-NW_f$       |        $-NW_b$        |   $-NW_g$   |        $-NW$         |

期末純資産は、期首純資産から現金の増減までの合計に等しい

### 1.1.4. ストックとフローの関係の整合性を示す表
期首の量とキャピタルゲインと関連するフローを足すと期末の量になる

| 資産もしくは負債/純資産の種類  | 期首の量  |   キャピタルゲイン    |     関連するフロー      | 期末の値 |
| :----------------------------- | :-------: | :-------------------: | :---------------------: | :------: |
| 資本(名目)                     | $K_{-1}$  |  $+\Delta p k_{-1}$   | $+p (-\delta k_{-1}+i)$ |   $K$    |
| 資本(実質)                     | $k_{-1}$  |                       |   $-\delta k_{-1}+i$    |   $k$    |
| 家計の預金残高                 | $M_{h-1}$ |                       |      $+\Delta M_h$      |  $M_h$   |
| 企業の預金残高                 | $M_{f-1}$ |                       |      $+\Delta M_f$      |  $M_f$   |
| 銀行の預金                     | $M_{-1}$  |                       |       $+\Delta M$       |   $M$    |
| 家計の借入金                   | $L_{h-1}$ |                       |      $+\Delta L_h$      |  $L_h$   |
| 企業の借入金                   | $L_{f-1}$ |                       |      $+\Delta L_f$      |  $L_f$   |
| 貸付金                         | $L_{-1}$  |                       |       $+\Delta L$       |   $L$    |
| 家計が保有する企業株式(時価)   | $E_{h-1}$ | $+\Delta p_e e_{h-1}$ |    $+p_e\Delta e_h$     |  $E_h$   |
| 家計が保有する企業株式(発行量) | $e_{h-1}$ |                       |      $+\Delta e_h$      |  $e_h$   |
| 銀行が保有する企業株式(時価)   | $E_{b-1}$ | $+\Delta p_e e_{b-1}$ |    $+p_e\Delta e_b$     |  $E_b$   |
| 銀行が保有する企業株式(発行量) | $e_{b-1}$ |                       |      $+\Delta e_b$      |  $e_b$   |
| 企業の資本金                   | $E_{-1}$  |                       |     $+p_e\Delta e$      |   $E$    |
| 株式発行量                     | $e_{-1}$  |                       |       $+\Delta e$       |   $e$    |
| 現金                           | $H_{-1}$  |                       |       $+\Delta H$       |   $H$    |

### 1.1.5. 表で示されない恒等式
- $p=p_{-1}+\Delta p$
- $p_e=p_{e-1}+\Delta p_e$

## 1.2. モデルの式一覧(計算する順)
- $p=\frac{(1+\nu_1+\nu_2\frac{L_{f-1}}{C_{-1}+G_{-1}})(W_{-1}+T_{v-1}+T_{c-1})}{u^T\gamma_1 k_{-1}}$
  - 価格を安定させたがること(塩沢)を反映した価格設定式はないか？
  - マークアップ率は価格競争と投資資金回収率を主な引数とする関数で書かれるのでは？集計的SFCモデルでは投資資金回収率に比例する項と定数項の和が価格マークアップの水準になるのが現実的か？
- $T_i=\tau_1 (W_{-1}+P_{h-1}+S_{-1})$
- $T_a=\tau_2(E_{h-1}+M_{h-1}-L_{h-1})$
- $T_v=\tau_3(C_{-1}+I_{-1}+G_{-1})$
- $T_c=\tau_4(C_{-1}+G_{-1}+I_{-1}-W_{-1}-T_v)$
- $G=(1+\beta)G_{-1}$
- $g=\frac{G}{p}$
- $C=\alpha_1(W_{-1}-T_a-T_i-r_L L_{h-1}+P_{h-1}+S_{-1}) + \alpha_2(E_{h-1}+M_{h-1}-L_{h-1})$
- $c=\frac{C}{p}$
- $A=A_{-1}(1 + \mu_1 + \mu_2 \frac{i_{-1}}{k_{-1}})$
  - 技術水準上昇率が投資の一次関数として決まる。
- $u = \frac{c+i+g}{\gamma_1 k_{-1}}$
- $i=\delta k_{-1} + (u_{-1}-u^T)\gamma_2 k_{-1} + \gamma_3\frac{M_{f-1}-L_{f-1}}{p}$
- $I=p i$
- $k=(1-\delta)k_{-1}+i$
- $K=p k$
- $v=u \frac{k_{-1}}{A}$
  - 労働稼働率は、技術水準(と労働者数の積)/消費および投資の実質需要　で定義する。
- $w=(1-\zeta_1)w_{-1} + \zeta_1 w_{-1}\{1 + \zeta_2(v - v^T)\}$
- $W=w(c_{-1}+g_{-1}+i_{-1})$
  - 家計側の賃金交渉力が高いときは名目賃金が上昇しやすく、逆に家計側の賃金交渉力が低いときは名目賃金が上がりにくい。家計側の賃金交渉力の相対的な強さは、労働稼働率の一次関数で決まることにする。
- $P=C+G+I-W-T_c-T_v-r_L L_{f-1}$
- $P_h=\max\{0, \theta_1(P-I)+\theta_2(M_{f-1}-L_{f-1})\}\frac{e_{h-1}}{e_{-1}}$
- $P_b=\max\{0, \theta_1(P-I)+\theta_2(M_{f-1}-L_{f-1})\}\frac{e_{b-1}}{e_{-1}}$
- $P_f=P-P_h-P_b$
- $S=\theta_3(r_L L_{-1}+P_b)+\theta_4 E_{b-1}$
- $NL_h=-C+W-T_i-T_a-r_L L_{h-1}+P_h+S$
- $NL_f=-I+P_f$
- $NL_b=r_L L_{-1}+P_b-S$
- $NL_g=-G+T_i+T_a+T_v+T_c$
- $L_h=\epsilon_1 NW_{h-1}+\epsilon_2 C$
- $\Delta L_h=L_h-L_{h-1}$
- $\Delta L_f=\max\{-L_{f-1}, (\lambda_3 + \lambda_4(\frac{P - P_f}{E_{h-1} + E_{b-1}} - r_L))(I+W+T_v+T_c+r_L L_{f-1} - \phi M_{f-1})\}$
  - もっとリアルな貸付金水準の行動方程式はないか？要調査
- $L_f=L_{f-1}+\Delta Lf$
- $L=L_h+L_f$
- $\Delta L=L-L_{-1}$
- $\Delta e=\frac{1}{p_{-1}}(1-\lambda_3 - \lambda_4(\frac{P - P_f}{E_{h-1} + E_{b-1}} - r_L))(I+W+T_v+T_c+r_L L_{f-1} - \phi M_{f-1})$
- $E=E_{-1}+p_{e-1}\Delta e$
- $e=e_{-1}+\Delta e$
- $E_h=(\lambda_1+\lambda_2\frac{P_{h-1}}{E_{h-1}})(NL_h+\Delta L_h+E_{h-1}+M_{h-1})$
- $E_b=(1-\lambda_5 r_L+\lambda_5 \frac{P_{b-1}}{E_{b-1}})(NL_b+E_{b-1})$
- $p_e=\frac{E_h+E_b}{e}$
- $e_h=\frac{E_h}{p_e}$
- $e_b=\frac{E_b}{p_e}$
- $\Delta e_h=e_h-e_{h-1}$
- $\Delta e_b=e_h-e_{b-1}$
- $\Delta M_h=NL_h-p_{e-1}\Delta e_h+\Delta L_h$
- $M_h=M_{h-1}+\Delta M_h$
- $\Delta M_f=NL_f+\Delta L_f+p_{e-1}\Delta e$
- $M_f=M_{f-1}+\Delta M_f$
- $M=M_h+M_f$
- $\Delta M=M-M_{-1}$
- $\Delta H = NL_b - p_{e-1}\Delta e_b + \Delta M - \Delta L$
- $H = H_{-1} + \Delta H$
- $NW_h=M_h-L_h+E_h$
- $NW_f=K+M_f-L_f-E$
- $NW_b=-M+L+E_b+H$
- $NW_g=H$
- $D=-E+E_h+E_b$
- $NW=NW_h+NW_f+NW_b+NW_g$

隠された恒等式として

- $K+D-NW=0$
- $NL_h+NL_f+NL_b+NL_g=0$
- $\Delta L=\Delta L_h+\Delta L_f$
- $\Delta M=\Delta M_h+\Delta M_f$
- $NL_g+\Delta H=0$

を使ってモデルの会計的一貫性を確認する

# 2. エージェント・ベース・モデル

AB_simulatio,n.jlで実装する予定

## 2.1. エージェント番号
|      | インデックス | 総数 |
| :--- | :----------: | :--: |
| 家計 |     $l$      | $L$  |
| 企業 |     $o$      | $O$  |
| 銀行 |     $n$      | $N$  |

変数の右上につけて、どのエージェントとどのエージェントの取引もしくは貸借関係なのかを区別する

## 2.2. 部門集計レベルの会計恒等式群
### 2.2.1. 取引フロー表
|                |                    家計                    |            企業(経常)             |            企業(資本)            |                    銀行                    |       統合政府       |  合計   |
| :------------- | :----------------------------------------: | :-------------------------------: | :------------------------------: | :----------------------------------------: | :------------------: | :-----: |
| 消費           |          $-\sum_l\sum_o C^{o,l}$           |      $+\sum_l\sum_o C^{o,l}$      |                                  |                                            |                      |   $0$   |
| 政府支出       |                                            |           $+\sum_o G^o$           |                                  |                                            |    $-\sum_o G^o$     |   $0$   |
| 投資           |                                            |           $+\sum_o I^o$           |          $-\sum_o I^o$           |                                            |                      |   $0$   |
| 賃金           |          $+\sum_l\sum_o W^{l,o}$           |      $-\sum_l\sum_o W^{l,o}$      |                                  |                                            |                      |   $0$   |
| 所得税         |              $-\sum_l T_i^l$               |                                   |                                  |                                            |   $+\sum_l T_i^l$    |   $0$   |
| 資産税         |              $-\sum_l T_a^l$               |                                   |                                  |                                            |   $+\sum_l T_a^l$    |   $0$   |
| 付加価値税     |                                            |          $-\sum_o T_v^o$          |                                  |                                            |   $+\sum_o T_v^o$    |   $0$   |
| 法人税         |                                            |          $-\sum_o T_c^o$          |                                  |                                            |   $+\sum_o T_c^o$    |   $0$   |
| 借入金利払     |     $-r_L \sum_l\sum_n L_{h-1}^{l,n}$      | $-r_L \sum_o\sum_n L_{f-1}^{o,n}$ |                                  |           $+r_L \sum_n L_{-1}^n$           |                      |   $0$   |
| 企業利潤       |         $+\sum_l\sum_o P_h^{l,o}$          |           $-\sum_o P^o$           |         $+\sum_o P_f^o$          |         $+\sum_n\sum_o P_b^{n,o}$          |                      |   $0$   |
| 銀行配当       |          $+\sum_l\sum_n S^{l,n}$           |                                   |                                  |          $-\sum_l\sum_n S^{l,n}$           |                      |   $0$   |
| [ 純貸出 ]     |            [ $\sum_l NL_h^l$ ]             |                                   |       [ $\sum_o NL_f^o$ ]        |            [ $\sum_n NL_b^n$ ]             |      [ $NL_g$ ]      | [ $0$ ] |
| 企業株式純発行 | $-\sum_l\sum_o p_{e-1}^o \Delta e_h^{o,l}$ |                                   |  $+\sum_o p_{e-1}^o \Delta e^o$  | $-\sum_n\sum_o p_{e-1}^o \Delta e_b^{o,n}$ |                      |   $0$   |
| 銀行株式純発行 |  $-\sum_l\sum_n p_{f-1}^n \Delta f^{n,l}$  |                                   |                                  |  $+\sum_l\sum_n p_{f-1}^n \Delta f^{n,l}$  |                      |   $0$   |
| 預金の移動     |      $-\sum_l\sum_n \Delta M_h^{n,l}$      |                                   | $-\sum_o\sum_n \Delta M_f^{n,o}$ |            $+\sum_n \Delta M^n$            |                      |   $0$   |
| 貸付金の移動   |      $+\sum_l\sum_n \Delta L_h^{l,n}$      |                                   | $+\sum_o\sum_n \Delta L_f^{o,n}$ |            $-\sum_n \Delta L^n$            |                      |   $0$   |
| 現金の移動     |                                            |                                   |                                  |            $-\sum_n \Delta H^n$            | $+\sum_n \Delta H^n$ |   $0$   |
| 合計           |                    $0$                     |                $0$                |               $0$                |                    $0$                     |         $0$          |         |

### 2.2.2. バランスシート表
|          |           家計            |           企業            |           銀行            |   統合政府    |     合計      |
| :------- | :-----------------------: | :-----------------------: | :-----------------------: | :-----------: | :-----------: |
| 資本     |                           |       $+\sum_o K^o$       |                           |               | $+\sum_o K^o$ |
| 預金     | $+\sum_l\sum_n M_h^{n,l}$ | $+\sum_o\sum_n M_f^{n,o}$ |       $-\sum_n M^n$       |               |      $0$      |
| 貸付金   | $-\sum_l\sum_n L_h^{l,n}$ | $-\sum_o\sum_n L_f^{o,n}$ |       $+\sum_n L^n$       |               |      $0$      |
| 企業株式 | $+\sum_l\sum_o E_h^{o,l}$ |       $-\sum_o E^o$       | $+\sum_n\sum_o E_b^{o,n}$ |               |    $+D_E$     |
| 銀行株式 | $+\sum_l\sum_o F_h^{o,l}$ |                           |       $-\sum_n F^n$       |               |    $+D_F$     |
| 現金     |                           |                           |       $+\sum_n H^n$       | $-\sum_n H^n$ |      $0$      |
| 純資産   |     $-\sum_l NW_h^l$      |     $-\sum_o NW_f^o$      |     $-\sum_n NW_b^n$      |    $-NW_g$    |     $-NW$     |
| 合計     |            $0$            |            $0$            |            $0$            |      $0$      |      $0$      |

### 2.2.3. 負債及び資産の増減と、純資産の増減の、整合性を表す表(Full-integration matrix)
|                          |                   家計                   |               企業               |                   銀行                   |       統合政府       |                  合計                   |
| :----------------------- | :--------------------------------------: | :------------------------------: | :--------------------------------------: | :------------------: | :-------------------------------------: |
| 期首純資産               |           $\sum_l NW_{h-1}^l$            |       $\sum_o NW_{f-1}^o$        |           $\sum_n NW_{b-1}^n$            |      $NW_{g-1}$      |                $NW_{-1}$                |
| 資本のキャピタルゲイン   |                                          |  $+\sum_o \Delta p^o k_{-1}^o$   |                                          |                      |      $+\sum_o \Delta p^o k_{-1}^o$      |
| 資本増減                 |                                          |     $-\sum_o p^o \Delta k^o$     |                                          |                      |        $-\sum_o p^o \Delta k^o$         |
| 預金の増減               |     $+\sum_l\sum_n \Delta M_h^{n,l}$     | $+\sum_o\sum_n \Delta M_f^{n,o}$ |           $-\sum_n \Delta M^n$           |                      |                   $0$                   |
| 貸付金の増減             |     $-\sum_l\sum_n \Delta L_h^{l,n}$     | $-\sum_o\sum_n \Delta L_h^{o,n}$ |           $+\sum_n \Delta L^n$           |                      |                   $0$                   |
| 企業株式キャピタルゲイン | $+\sum_l\sum_o \Delta p_e e_{h-1}^{o,l}$ |                                  | $+\sum_n\sum_o \Delta p_e e_{b-1}^{o,n}$ |                      |      $+\sum_n \Delta p_e e_{-1}^n$      |
| 企業株式純発行           |   $+\sum_l\sum_o p_e \Delta e_h^{o,l}$   |     $+\sum_o p_e \Delta e^o$     |   $+\sum_n\sum_o p_e \Delta e_b^{o,n}$   |                      |                   $0$                   |
| 銀行株式キャピタルゲイン | $+\sum_l\sum_n \Delta p_f f_{-1}^{n,l}$  |                                  |                                          |                      | $+\sum_l\sum_n \Delta p_f f_{-1}^{n,l}$ |
| 銀行株式純発行           |    $+\sum_l\sum_n p_f \Delta f^{n,l}$    |                                  |    $-\sum_l\sum_n p_f \Delta f^{n,l}$    |                      |                                         |
| 現金の増減               |                                          |                                  |           $+\sum_n \Delta H^n$           | $-\sum_n \Delta H^n$ |                   $0$                   |
| 純資産                   |             $-\sum_l NW_h^l$             |         $-\sum_o NW_f^o$         |             $-\sum_n NW_b^n$             |       $-NW_g$        |                  $-NW$                  |

期末純資産は、期首純資産から現金の増減までの合計に等しい

## 2.3. エージェントレベルの会計恒等式群
経済主体にとっての借方と貸方のバランスに相当する

### 2.3.1. フローの一貫性(貸借対照表及びキャッシュ・フロー計算書のようなもの)の表
|                |                 家計                 |         企業(経常)          |         企業(資本)         |                 銀行                 |       統合政府       |
| :------------- | :----------------------------------: | :-------------------------: | :------------------------: | :----------------------------------: | :------------------: |
| 消費           |          $-\sum_o C^{o,l}$           |      $+\sum_l C^{o,l}$      |                            |                                      |                      |
| 政府支出       |                                      |           $+G^o$            |                            |                                      |    $-\sum_o G^o$     |
| 投資           |                                      |           $+I^o$            |           $-I^o$           |                                      |                      |
| 賃金           |          $+\sum_o W^{l,o}$           |      $-\sum_l W^{l,o}$      |                            |                                      |                      |
| 所得税         |               $-T_i^l$               |                             |                            |                                      |   $+\sum_l T_i^l$    |
| 資産税         |               $-T_a^l$               |                             |                            |                                      |   $+\sum_l T_a^l$    |
| 付加価値税     |                                      |          $-T_v^o$           |                            |                                      |   $+\sum_o T_v^o$    |
| 法人税         |                                      |          $-T_c^o$           |                            |                                      |   $+\sum_o T_c^o$    |
| 借入金利払     |     $-r_L \sum_n L_{h-1}^{l,n}$      | $-\sum_n r_L L_{f-1}^{o,n}$ |                            |           $+r_L L_{-1}^n$            |                      |
| 企業利潤       |         $+\sum_o P_h^{l,o}$          |           $-P^o$            |          $+P_f^o$          |         $+\sum_o P_b^{n,o}$          |                      |
| 銀行配当       |          $+\sum_n S^{l,n}$           |                             |                            |          $-\sum_l S^{l,n}$           |                      |
| [ 純貸出 ]     |             [ $NL_h^l$ ]             |                             |        [ $NL_f^o$ ]        |             [ $NL_b^n$ ]             |      [ $NL_g$ ]      |
| 企業株式純発行 | $-\sum_o p_{e-1}^o \Delta e_h^{o,l}$ |                             |  $+p_{e-1}^o \Delta e^o$   | $-\sum_o p_{e-1}^o \Delta e_b^{o,n}$ |                      |
| 銀行株式純発行 |  $-\sum_n p_{f-1}^n \Delta f^{n,l}$  |                             |                            |  $+\sum_l p_{f-1}^n \Delta f^{n,l}$  |                      |
| 預金の移動     |      $-\sum_n \Delta M_h^{n,l}$      |                             | $-\sum_n \Delta M_f^{n,o}$ |            $+\Delta M^n$             |                      |
| 貸付金の移動   |      $+\sum_n \Delta L_h^{l,n}$      |                             | $+\sum_n \Delta L_f^{o,n}$ |            $-\Delta L^n$             |                      |
| 現金の移動     |                                      |                             |                            |            $-\Delta H^n$             | $+\sum_n \Delta H^n$ |
| 合計           |                 $0$                  |             $0$             |            $0$             |                 $0$                  |         $0$          |

### 2.3.2. ストックの一貫性(バランスシート)の表
|          |        家計         |        企業         |        銀行         |   統合政府    |
| :------- | :-----------------: | :-----------------: | :-----------------: | :-----------: |
| 資本     |                     |       $+K^o$        |                     |               |
| 預金     | $+\sum_n M_h^{n,l}$ | $+\sum_n M_f^{n,o}$ |       $-M^n$        |               |
| 貸付金   | $-\sum_n L_h^{l,n}$ | $-\sum_n L_f^{o,n}$ |       $+L^n$        |               |
| 企業株式 | $+\sum_o E_h^{o,l}$ |       $-E^o$        | $+\sum_o E_b^{o,n}$ |               |
| 銀行株式 |  $+\sum_o F^{n,l}$  |                     |       $-F^n$        |               |
| 現金     |                     |                     |       $+H^n$        | $-\sum_n H^n$ |
| 純資産   |      $-NW_h^l$      |      $-NW_f^o$      |      $-NW_b^n$      |    $-NW_g$    |
| 合計     |         $0$         |         $0$         |         $0$         |      $0$      |

### 2.3.3. 負債及び資産の増減と、純資産の増減の、整合性を表す表(Full-integration matrix)
|                          |                家計                |            企業            |                銀行                |       統合政府       |
| :----------------------- | :--------------------------------: | :------------------------: | :--------------------------------: | :------------------: |
| 期首純資産               |            $NW_{h-1}^l$            |        $NW_{f-1}^o$        |            $NW_{b-1}^n$            |      $NW_{g-1}$      |
| 資本のキャピタルゲイン   |                                    |   $+\Delta p^o k_{-1}^o$   |                                    |                      |
| 資本増減                 |                                    |     $+p^o \Delta k^o$      |                                    |                      |
| 預金の増減               |     $+\sum_n \Delta M_h^{n,l}$     | $+\sum_n \Delta M_f^{n,o}$ |           $-\Delta M^n$            |                      |
| 貸付金                   |     $-\sum_n \Delta L_h^{l,n}$     | $-\sum_n \Delta L_h^{o,n}$ |           $+\Delta L^n$            |                      |
| 企業株式キャピタルゲイン | $+\sum_o \Delta p_e e_{h-1}^{o,l}$ |                            | $+\sum_o \Delta p_e e_{b-1}^{o,n}$ |                      |
| 企業株式純発行           |   $+\sum_o p_e \Delta e_h^{o,l}$   |     $-p_e \Delta e^o$      |   $+\sum_o p_e \Delta e_b^{o,n}$   |                      |
| 銀行株式キャピタルゲイン | $+\sum_n \Delta p_f f_{-1}^{n,l}$  |                            |                                    |                      |
| 銀行株式純発行           |    $+\sum_n p_f \Delta f^{n,l}$    |                            |    $-\sum_l p_f \Delta f^{n,l}$    |                      |
| 現金の増減               |                                    |                            |        $+\sum_n \Delta H^n$        | $-\sum_n \Delta H^n$ |
| 期末純資産               |             $-NW_h^l$              |         $-NW_f^o$          |             $-NW_b^n$              |       $-NW_g$        |

期末純資産は、期首純資産から現金の増減までの合計に等しい

### 2.3.4. ストックとフローの関係の整合性を示す表
期首の量とキャピタルゲインと関連するフローを足すと期末の量になる

| 資産もしくは負債/純資産の種類  |    期首の量     |      キャピタルゲイン       |        関連するフロー         |  期末の値   |
| :----------------------------- | :-------------: | :-------------------------: | :---------------------------: | :---------: |
| 資本(名目)                     |   $K_{-1}^o$    |   $+\Delta p^o k_{-1}^o$    | $+p^o (-\delta k_{-1}^o+i^o)$ |    $K^o$    |
| 資本(実質)                     |   $k_{-1}^o$    |                             |    $-\delta k_{-1}^o+i^o$     |    $k^o$    |
| 家計の預金残高                 | $M_{h-1}^{n,l}$ |                             |      $+\Delta M_h^{n,l}$      | $M_h^{n,l}$ |
| 企業の預金残高                 | $M_{f-1}^{n,o}$ |                             |      $+\Delta M_f^{n,o}$      | $M_f^{n,o}$ |
| 銀行の預金                     |   $M_{-1}^n$    |                             |         $+\Delta M^n$         |    $M^n$    |
| 家計の借入金                   | $L_{h-1}^{n,l}$ |                             |      $+\Delta L_h^{n,l}$      | $L_h^{n,l}$ |
| 企業の借入金                   | $L_{f-1}^{n,o}$ |                             |      $+\Delta L_f^{n,o}$      | $L_f^{n,o}$ |
| 貸付金                         |   $L_{-1}^n$    |                             |         $+\Delta L^n$         |    $L^n$    |
| 家計が保有する企業株式(時価)   | $E_{h-1}^{o,l}$ | $+\Delta p_e e_{h-1}^{o,l}$ |    $+p_e\Delta e_h^{o,l}$     | $E_h^{o,l}$ |
| 家計が保有する企業株式(発行量) | $e_{h-1}^{o,l}$ |                             |      $+\Delta e_h^{o,l}$      | $e_h^{o,l}$ |
| 銀行が保有する企業株式(時価)   | $E_{b-1}^{o,n}$ | $+\Delta p_e e_{b-1}^{o,n}$ |    $+p_e\Delta e_b^{o,n}$     | $E_b^{o,n}$ |
| 銀行が保有する企業株式(発行量) | $e_{b-1}^{o,n}$ |                             |      $+\Delta e_b^{o,n}$      | $e_b^{o,n}$ |
| 企業の資本金                   |   $E_{-1}^o$    |                             |       $+p_e\Delta e^o$        |    $E^o$    |
| 企業株式発行量                 |   $e_{-1}^o$    |                             |         $+\Delta e^o$         |    $e^o$    |
| 家計が保有する銀行株式(時価)   | $F_{h-1}^{n,l}$ | $+\Delta p_f f_{h-1}^{n,l}$ |    $+p_f\Delta f_h^{n,l}$     | $F_h^{n,l}$ |
| 家計が保有する銀行株式(発行量) | $f_{-1}^{n,l}$  |                             |       $+\Delta f^{n,l}$       |  $f^{n,l}$  |
| 銀行の資本金                   |   $F_{-1}^n$    |                             |    $+p_f\Delta f_h^{n,l}$     |    $F^n$    |
| 現金                           |   $H_{-1}^n$    |                             |         $+\Delta H^n$         |    $H^n$    |

表の内容＋αを式にするとこうなる
- $K^o=K_{-1}^o+\Delta p^o k_{-1}^o+p^o \Delta k^o=K_{-1}^o+\Delta p^o k_{-1}^o+p^o (-\delta k_{-1}^o+i^o)$
- $k^o=k_{-1}^o+\Delta k^o=k_{-1}^o+p^o (-\delta k_{-1}^o+i^o)$
- $M_h^{n,l}=M_{h-1}^{n,l}+\Delta M_h^{n,l}$
- $M_f^{n,o}=M_{f-1}^{n,o}+\Delta M_f^{n,o}$
- $M^n=M_{-1}^n+\Delta M^n$
- $L_h^{n,l}=L_{h-1}^{n,l}+\Delta L_h^{n,l}$
- $L_f^{n,o}=L_{f-1}^{n,o}+\Delta L_f^{n,o}$
- $L^n=L_{-1}^n+\Delta L^n$
- $E_h^{o,l}=E_{h-1}^{o,l}+\Delta p_e e_{h-1}^{o,l}+p_e\Delta e_h^{o,l}$
- $e_h^{o,l}=e_{h-1}^{o,l}+\Delta e_h^{o,l}$
- $E_b^{o,n}=E_{b-1}^{o,n}+\Delta p_e e_{b-1}^{o,n}+p_e\Delta e_b^{o,n}$
- $e_b^{o,n}=e_{b-1}^{o,n}+\Delta e_b^{o,n}$
- $E^o=E_{-1}^o+p_e\Delta e^o$
- $e^o=e_{-1}^o+\Delta e^o$
- $F_h^{n,l}=F_{h-1}^{n,l}+\Delta p_f f_{h-1}^{n,l}+p_f\Delta f_h^{n,l}$
- $f^{n,l}=f_{-1}^{n,l}+\Delta f^{n,l}$
- $F^n=F_{-1}^n+p_f\Delta f_h^{n,l}$
- $H^n=H_{-1}^n+\Delta H^n$

## 2.4. 表で示されない恒等式
- $\sum_l P_h^{o,l}-P^o+P_f^o+\sum_n P_b^{o,n}=0$
- $e^o=\sum_l e_h^{o,l} + \sum_n e_b^{o,n}$
- $-\sum_l \Delta e_h^{o,l}+\Delta e^o-\sum_n \Delta e_b^{o,n}=0$
- $-\sum_l M_h^{n,l}-\sum_o M_f^{n,o}+M^n=0$
- $-\sum_l \Delta M_h^{n,l}-\sum_o \Delta M_f^{n,o}+\Delta M^n=0$
- $+\sum_l L_h^{n,l}+\sum_o L_f^{n,o}-L^n=0$
- $+\sum_l \Delta L_h^{n,l}+\sum_o \Delta L_f^{n,o}-\Delta L^n=0$
- $p^o=p_{-1}^o+\Delta p^o$
- $p_e^o=p_{e-1}^o+\Delta p_e^o$
- $p_f^n=p_{f-1}^n+\Delta p_f^n$




## 2.5. モデルの式とアルゴリズム
- $p^o=\lambda_p \frac{(1+\nu_1+\nu_2\frac{\sum_n L_{f-1}^{o,n}}{\sum_l C_{-1}^{o,l}+ G_{-1}^o})(\sum_l W_{-1}^{l,o}+T_{v-1}^o+T_{c-1}^o+\delta k_{-1}^o)}{u^T\gamma_1 k_{-1}^o} + (1-\lambda_p)\nu_3(p_{-1}^o-\overline{p}_{-1})$
  - マークアップ率は価格競争と投資資金回収率を主な引数とする関数で書かれるのでは？
  - 価格競争の成分は？ $\overline{p}_{-1}=\frac{\sum_o(\sum_l C^{o,l}_{-1}+G_{-1}^o)}{\sum_o(\sum_l c^{o,l}_{-1}+g_{-1}^o)}$ に近づくように $\nu_3(p_{-1}^o-\overline{p}_{-1})$ の項を追加する？
  - 収穫逓増の効果を入れるべき？入れるとしたらどうやって？
  - $\lambda_p$を適応ラメータにして、消費者が価格差に対しする商品の乗り換え速度に適応するようにしたほうが良いかも
- $T_i^l=\tau_1 \sum_l (W_{-1}^{l,o}+P_{h-1}^{l,o}+S_{-1}^{l,n})$
- $T_a^l=\tau_2(\sum_o E_{h-1}^{o,l}+\sum_n M_{h-1}^{n,l}-\sum_n L_{h-1}^{l,n})$
- $T_v^o=\tau_3(\sum_l C_{-1}^{o,l}+I_{-1}^o+G_{-1}^o)$
- $T_c^o=\tau_4(\sum_l C_{-1}^{o,l}+G_{-1}^o+I_{-1}^o-\sum_l W_{-1}^{l,o}-T_v^o)$
- 政府支出のアルゴリズム　govExpAlg2.jl のアルゴリズムを採用
  - $G_{sum}=\sum G^o=G_0(1+\beta_1)^{(t-1)}$
  - $x^o=\max[1, x_{-1}^o \{1+\beta_2(-1+\beta_3 randn)\}]$
  - $z^o=\max(0, x^o-\beta_4)$
  - $G^o=\frac{z^o G_{sum}}{\sum z^o}$
  - 下限付きべき分布みたいなものを作って(xがこれにあたる)、一定の値を差し引いた後(zがこれに当たる)に、政府支出総額に合わせて大きくする。
- $g^o=\frac{G^o}{p^o}$
- $\sum_o C^{o,l}=\alpha_1(\sum_o W_{-1}^{o,l}-T_a^l-T_i^l-r_L \sum_n L_{h-1}^{l,n}+\sum_o P_{h-1}^{o,l}+\sum_n S_{-1}^{n,l}) + \alpha_2(\sum_o E_{h-1}^{o,l}+\sum_n M_{h-1}^{n,l}-\sum_n L_{h-1}^{n,l})$
- 消費財の購入先のアルゴリズム
  - 前回消費財を買った企業$o$から別の企業に乗り換える確率は
    - $rand() < \alpha_3 \frac{p^{o}-\overline{p}}{\overline{p}}$
    - とする。乗り換え先は「前期の消費財生産量($\sum_l c_{-1}^{o,l}$)＋前期の総消費量の一定割合($\alpha_4 \sum_o \sum_l c_{-1}^{o,l}$)」に比例する確率でランダムに決める
  - 今期の消費財購入先の企業を$o'$と書くことにすると、今期の消費水準は、
    - $c^{o'l}=\frac{\sum_o C^{o,l}}{p^{o'}}$
    - で決める。
    - $o'$以外のすべての企業の商品は消費しない($c^{o,l}=0 \ \ \ \{ o \mid o \neq o'\}$)
- $C^{o,l}=p^o c^{o,l}$
- $A^o=A_{-1}^o(1 + \mu_1 + \mu_2 \frac{i_{-1}^o}{k_{-1}^o})$
  - 技術水準上昇率が投資の一次関数として決まる。
- $u^o = \frac{\sum_l c^{o,l}+i^o+g^o}{\gamma_1 k_{-1}^o}$
- $i^o=\delta k_{-1}^o + (u_{-1}^o-u^T)\gamma_2 k_{-1}^o + \gamma_3\frac{\sum_n M_{f-1}^{n,o}-\sum_n L_{f-1}^{o,n}}{p^o}$
- $I^o=p^o i^o$
- $k^o=(1-\delta)k_{-1}^o+i^o$
- $K^o=p^o k^o$
- 就業・失業判定アルゴリズム　外部にファイルを作る
  - if $\sum_o EMP_{-1}^{o,l}=1$
    - if $rand() < \zeta_3$
      - $EMP^{o,l}=0 \ \ \  \forall o$
      - $\sum_o EMP^{o,l}=0$
      - 失業者への社会保障をモデル含めることを検討
      - $w^l=(1-\zeta_1)w_{-1}^l + \zeta_1 w_{-1}[1 - \zeta_3 abs\{rand()\}]$
    - else
      - $EMP^{o,l}=EMP_{-1}^{o,l} \ \ \  \forall o$
      - $\sum_o EMP^{o,l}=1$
      - $w^l=(1-\zeta_1)w_{-1}^l + \zeta_1 w_{-1}[1 + \zeta_2 abs\{rand()\}]$
  - else
    - 企業が Int64($\max[0, \frac{1}{A^o}\{u_{-1}^o k_{-1}^o-A^o \sum_l (w_{-1}^{l,o}>0)\}]$) 人の求人を出す
    - 求人数に比例する確率で失業者は応募する。応募の中から定員までランダムに雇用する
    - $o'$で雇用が決まった場合
      - $w^l=(1-\zeta_1)w_{-1}^l + \zeta_1 w_{-1}[1 + \zeta_2 abs\{rand()\}]$
      - $EMP^{o'l}=1$
      - $EMP^{o,l}=0 \ \ \ (o\neq o')$
    - 失業が続く場合
      - $w^l=(1-\zeta_1)w_{-1}^l + \zeta_1 w_{-1}[1 - \zeta_3 abs\{rand()\}]$
- $W^{l,o}=w^l EMP^{o,l}$
- $P^o=\sum_l C^{o,l}+G^o+I^o-\sum_l W^{l,o}-T_c^o-T_v^o-r_L \sum_n L_{f-1}^{o,n}$
- $P_h^{l,o}=\max\{0, \theta_1(P^o-I^o)+\theta_2(\sum_n M_{f-1}^{o,n}-\sum_n L_{f-1}^{o,n})\}\frac{e_{h-1}^{o,l}}{e_{-1}^o}$
- $P_b^{n,o}=\max\{0, \theta_1(P^o-I^o)+\theta_2(\sum_n M_{f-1}^{o,n}-\sum_n L_{f-1}^{o,n})\}\frac{e_{b-1}^{o,n}}{e_{-1}^o}$
- $P_f^o=P^o-\sum_l P_h^{l,o}-\sum_n P_b^{n,o}$
- $S^{l,n}=\{\theta_3(r_L L_{-1}^n+\sum_o P_b^{n,o})+\theta_4 \sum_o E_{b-1}^{o,n}\}\frac{f_{h-1}^{l,n}}{f_{-1}^n}$
- $NL_h^l=-\sum_o C^{o,l}+\sum_o W^{l,o}-T_i^l-T_a^l-r_L L_{h-1}^{l,n}+\sum_o P_h^{l,o}+\sum_n S^{l,n}$
- $NL_f^o=-I^o+P_f^o$
- $NL_b^n=r_L L_{-1}^n+\sum_o P_b^{o,n}-\sum_l S^{l,n}$
- $NL_g=-\sum_o G^o+\sum_l T_i^l+\sum_l T_a^l+\sum_o T_v^o+\sum_o T_c^o$
- $\sum_o L_h^{l,n}=\epsilon_1 NL_h^l+\epsilon_2 \sum_o C^{o,l}$
  - 家計は借入先をどう選ぶ？
    - 既存の借入先があればそれを継続する
    - 既存の借入先がなければ、$NW_b^n$に比例する確率で借入先に選ぶ
    - 本当は銀行ごとに信用スコアを計算して金利が安いところから借り入れようとするとか、既存の借入先を維持する傾向を持つとか、借入先を2つ以上にすることもあるとか、そういう効果を入れたいが、モデルの複雑さを抑えるためのアドホックな仮定として導入することにする。
- $L_h^{l,n}=$
  - ここのアルゴリズムも別でファイルに書く。
- $\Delta L_h^{l,n}=L_h^{l,n}-L_{h-1}^{l,n}$
- $\sum_n \Delta L_f^{o,n}=\max\{-\sum_n L_{f-1}^{o,n}, (\lambda_3 + \lambda_4(\frac{P^o - P_f^o}{\sum_l E_{h-1}^{o,l} + \sum_l E_{b-1}^{o,l}} - r_L))(I^o+\sum_l W^{l,o}+T_v^o+T_c^o+r_L \sum_n L_{f-1}^{l,n} - \phi \sum_n M_{f-1}^{n,l})\}$
  - もっとリアルな貸付金水準の行動方程式はないか？要調査
- $L_f^{o,n}=$
  - $L_h^{l,n}$ と同じ方法で振り分ける
- 
- ここからAB化の作業を再開する
- 
- $L_f=L_{f-1}+\Delta L_f$
- $L=L_h+L_f$
- $\Delta L=L-L_{-1}$
- $\Delta e=\frac{1}{p_{-1}}(1-\lambda_3 - \lambda_4(\frac{P - P_f}{E_{h-1} + E_{b-1}} - r_L))(I+W+T_v+T_c+r_L L_{f-1} - \phi M_{f-1})$
- $E=E_{-1}+p_{e-1}\Delta e$
- $e=e_{-1}+\Delta e$
- $E_h=(\lambda_1+\lambda_2\frac{P_{h-1}}{E_{h-1}})(NL_h+\Delta L_h+E_{h-1}+M_{h-1})$
- $E_b=(1-\lambda_5 r_L+\lambda_5 \frac{P_{b-1}}{E_{b-1}})(NL_b+E_{b-1})$
- $p_e=\frac{E_h+E_b}{e}$
- $e_h=\frac{E_h}{p_e}$
- $e_b=\frac{E_b}{p_e}$
- $\Delta e_h=e_h-e_{h-1}$
- $\Delta e_b=e_h-e_{b-1}$
- $\Delta M_h=NL_h-p_{e-1}\Delta e_h+\Delta L_h$
- $M_h=M_{h-1}+\Delta M_h$
  - インターバンク市場のことも考えなきゃ
- $\Delta M_f=NL_f+\Delta L_f+p_{e-1}\Delta e$
- $M_f=M_{f-1}+\Delta M_f$
- $M=M_h+M_f$
- $\Delta M=M-M_{-1}$
- $\Delta H = NL_b - p_{e-1}\Delta e_b + \Delta M - \Delta L$
- $H = H_{-1} + \Delta H$
- $NW_h=M_h-L_h+E_h$
- $NW_f=K+M_f-L_f-E$
- $NW_b=-M+L+E_b+H$
- $NW_g=H$
- $D=-E+E_h+E_b$
- $NW=NW_h+NW_f+NW_b+NW_g$
- 倒産と新規参入のアルゴリズムも追加する必要がある

隠された恒等式として

- $K+D-NW=0$
- $NL_h+NL_f+NL_b+NL_g=0$
- $\Delta L=\Delta L_h+\Delta L_f$
- $\Delta M=\Delta M_h+\Delta M_f$
- $NL_g+\Delta H=0$

を使ってモデルの会計的一貫性を確認する

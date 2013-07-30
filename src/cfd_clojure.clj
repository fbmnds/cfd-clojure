(ns cfd-clojure
  (:use (incanter core charts)
        (clojure inspector)))

;; in O(n) on lazy-seqs
;;
(defn cum-fn
  ([s cfn]
     (cond (empty? s) nil
           (empty? (rest s)) (list (first s))
           :else (lazy-seq (cons (first s) (cum-fn (first s) (rest s) cfn)))))
  ([x s cfn]
     (cond (empty? s) nil
           (empty? (rest s)) (list (cfn x (first s)))
           :else (lazy-seq (cons (cfn x (first s)) (cum-fn (cfn x (first s)) (rest s) cfn))))))


(defn u_i [u_ni u_ni-1 c dx dt]
  (- u_ni (/ (* c dt (- u_ni u_ni-1)) dx)))

;; ~>   (reduce #(u_i %2 %1 c dx dt) u)

(defn linear-convection [u0 c nx dx nt dt]
  (let [u0_0 (first u0)]
    (loop [u u0
           n 0]
      (if (= n nt)
        u
        (recur (conj (map #(u_i (second %) (first %) c dx dt) (partition 2 1 u)) u0_0)
               (inc n))))))

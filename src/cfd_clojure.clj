(ns cfd-clojure
  (:use (incanter core charts)))



;;(defn u_i [u_ni u_ni-1 c dx dt]
;;  (- u_ni (/ (* c dt (- u_ni u_ni-1)) dx)))
;;
;;
;;
;;(defn linear-convection [u0 c nx dx nt dt]
;;  (let [u0_0 (first u0)]
;;    (loop [u u0
;;           n 0]
;;      (if (= n nt)
;;        u
;;        (recur (conj (map #(u_i (second %) (first %) c dx dt) (partition 2 1 u)) u0_0)
;;               (inc n))))))



(defn linear-convection [m [u_ni-1 u_ni]]
  (- u_ni (/ (* (:c m) (:dt m) (- u_ni u_ni-1)) (:dx m))))



;;(defn discretize
;;  ([f m u]
;;     (let [r (conj u (conj (map #(f m %) (partition {:x-steps m} 1 u)) (first u)))]
;;       (lazy-seq (cons (discretize f m (last r) r)))
;;  ([f m u c]
;;     (let [r (conj (map #(f m %) (partition {:x-steps m} 1 u)) (first u))]
;;       (lazy-seq (cons (discretize f m r c) c)))))))

(defn discretize [f m u]
  (let [g (fn [u] (conj (map #(f m %) (partition (:x-steps m) 1 u)) (first u)))]
    (iterate g u)))

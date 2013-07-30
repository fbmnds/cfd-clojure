(ns cfd-clojure
  (:use (incanter core charts)))



(defn u_i [u_ni u_ni-1 c dx dt]
  (- u_ni (/ (* c dt (- u_ni u_ni-1)) dx)))



(defn linear-convection [u0 c nx dx nt dt]
  (let [u0_0 (first u0)]
    (loop [u u0
           n 0]
      (if (= n nt)
        u
        (recur (conj (map #(u_i (second %) (first %) c dx dt) (partition 2 1 u)) u0_0)
               (inc n))))))

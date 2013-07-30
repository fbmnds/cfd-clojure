(defproject cfd-clojure "0.0.1-SNAPSHOT"
  :description "CFD Python port to Clojurex"
  :dependencies [[org.clojure/clojure "1.5.1"]

                 [incanter/incanter-core "1.5.0-SNAPSHOT"]
                 [incanter/incanter-io "1.5.0-SNAPSHOT"]
                 [incanter/incanter-charts "1.5.0-SNAPSHOT"]
                 [incanter/incanter-mongodb "1.5.0-SNAPSHOT"]
                 [incanter/incanter-pdf "1.5.0-SNAPSHOT"]
                 [incanter/incanter-latex "1.5.0-SNAPSHOT"]
                 [incanter/incanter-excel "1.5.0-SNAPSHOT"]
                 [incanter/incanter-sql "1.5.0-SNAPSHOT"]
                 [incanter/incanter-zoo "1.5.0-SNAPSHOT"]
                 [swingrepl "1.3.0"
                  :exclusions [org.clojure/clojure
                               org.clojure/clojure-contrib]]
                 [jline "0.9.94"]]
  :profiles {:dev {:dependencies [[midje "1.5.1"]]}}
  :repl-options {:port 4555})

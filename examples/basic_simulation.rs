use stochastic_gbm::gbm::GeometricBrownianMotion;

fn main() {
    // Daha gerçekçi piyasa parametreleri:
    // mu = 0.05 (5% beklenen getiri)
    // sigma = 0.3 (30% volatilite)
    // 10 paths (örnek yollar)
    // 252 steps (bir yıldaki işlem günü sayısı)
    let mut gbm = GeometricBrownianMotion::new(0.05, 0.3, 10, 252, 1.0, 100.0)
        .expect("Failed to create GBM simulator");
    
    let result = gbm.simulate()
        .expect("Failed to run simulation");
    
    // İlk iki path'i gösterelim
    for path_idx in 0..2 {
        println!("\nPath {}:", path_idx + 1);
        for (i, price) in result.paths[path_idx].iter().enumerate() {
            if i % 21 == 0 { // Her ay (yaklaşık 21 işlem günü)
                println!("t={:.2}: ${:.2}", i as f64 / 252.0, price);
            }
        }
    }

    // Son zaman noktası için istatistikleri gösterelim
    let (mean, std_dev, min, max) = GeometricBrownianMotion::get_statistics(&result, 252);
    println!("\nFinal Statistics (t=1.0):");
    println!("Mean Price: ${:.2}", mean);
    println!("Standard Deviation: ${:.2}", std_dev);
    println!("Min Price: ${:.2}", min);
    println!("Max Price: ${:.2}", max);
} 
"use strict"

document.addEventListener('DOMContentLoaded', function () {
        var ctx = document.getElementById('id_batteryChargeChart').getContext('2d');

        var myChart = new Chart(ctx, {
            type: 'line',
            data: {
                labels: P_ch_chart_labels,
                datasets: [{
                    label: 'Battery Charge Rate (W)',
                    data: P_ch_chart_data,
                    fill: false,
                    borderColor: 'rgb(75, 192, 223)',
                    tension: 0.1
                }]
            },
            options: {
                responsive: true,
                maintainAspectRatio: false,
                scales: {
                    y: {
                        beginAtZero: true
                    }
                }
            }
        });

        ctx = document.getElementById('id_batteryDischargeChart').getContext('2d');
        myChart = new Chart(ctx, {
            type: 'line',
            data: {
                labels: P_d_chart_labels,
                datasets: [{
                    label: 'Battery Discharge (W)',
                    data: P_d_chart_data,
                    fill: false,
                    borderColor: 'rgb(75, 192, 223)',
                    tension: 0.1
                }]
            },
            options: {
                responsive: true,
                maintainAspectRatio: false,
                scales: {
                    y: {
                        beginAtZero: true
                    }
                }
            }
        });

        ctx = document.getElementById('id_batteryStateofCharge').getContext('2d');
        myChart = new Chart(ctx, {
            type: 'line',
            data: {
                labels: B_chart_labels,
                datasets: [{
                    label: 'Battery State of Charge (pu)',
                    data: B_chart_data,
                    fill: false,
                    borderColor: 'rgb(75, 192, 223)',
                    tension: 0.1
                }]
            },
            options: {
                responsive: true,
                maintainAspectRatio: false,
                scales: {
                    y: {
                        beginAtZero: true
                    }
                }
            }
        });

        ctx = document.getElementById('id_slackGeneration').getContext('2d');
        myChart = new Chart(ctx, {
            type: 'line',
            data: {
                labels: sg_chart_labels,
                datasets: [{
                    label: 'Imported Power (W)',
                    data: sg_chart_data,
                    fill: false,
                    borderColor: 'rgb(75, 192, 223)',
                    tension: 0.1
                }]
            },
            options: {
                responsive: true,
                maintainAspectRatio: false,
                scales: {
                    y: {
                        beginAtZero: true
                    }
                }
            }
        });

        ctx = document.getElementById('id_renewGeneration').getContext('2d');
        myChart = new Chart(ctx, {
            type: 'line',
            data: {
                labels: renew_chart_labels,
                datasets: [{
                    label: 'Renewable Generation (W)',
                    data: renew_chart_data,
                    fill: false,
                    borderColor: 'rgb(75, 192, 223)',
                    tension: 0.1
                }]
            },
            options: {
                responsive: true,
                maintainAspectRatio: false,
                scales: {
                    y: {
                        beginAtZero: true
                    }
                }
            }
        });


    });